import io
import re
import torch
import dash
import plotly.express as px

from dash import dcc, html, Input, Output, State
from transformers import AutoModelForCausalLM, AutoTokenizer

# Global variables to store model and tokenizer
model = None
tokenizer = None

def generate_sequences(prompt, min_len, max_len, rep_pen, num_seq, top_k=950):
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model.to(device)
    inputs = tokenizer([prompt], return_tensors="pt").to(device)
    generated_ids = model.generate(
        **inputs,
        min_length=min_len,
        max_length=max_len,
        do_sample=True,
        top_k=top_k,
        repetition_penalty=rep_pen,
        num_return_sequences=num_seq,
        eos_token_id=2,
        truncation=True
    )

    # Decode and clean outputs
    generated_sequences = []
    for seq_ids in generated_ids:
        seq = tokenizer.decode(seq_ids, skip_special_tokens=True)
        seq = seq.replace(' ', '')
        generated_sequences.append(seq)
    return generated_sequences


def to_fasta(seqs: list[str]) -> str:
    lines = []
    for i, s in enumerate(seqs, 1):
        lines.append(f">sarita_seq_{i}")
        lines.extend(re.findall(".{1,70}", s))
    return "\n".join(lines) + "\n"

# Initialize Dash app
app = dash.Dash(__name__)

# Define styles
CARD_STYLE = {
    'backgroundColor': '#ffffff',
    'padding': '20px',
    'borderRadius': '10px',
    'boxShadow': '0 4px 6px rgba(0, 0, 0, 0.1)',
    'marginBottom': '20px'
}
CONTAINER_STYLE = {
    'maxWidth': '800px',
    'margin': '40px auto',
    'fontFamily': 'Arial, sans-serif',
    'color': '#333'
}
LABEL_STYLE = {'fontWeight': 'bold', 'marginTop': '10px'}

app.layout = html.Div(style=CONTAINER_STYLE, children=[
    html.H1("SARITA Spike Protein Generator", style={'textAlign': 'center', 'marginBottom': '30px'}),

    # Model info card
    html.Div(style=CARD_STYLE, children=[
        html.P("SARITA is a Large Language Model (LLM) specialized in generating synthetic S1 subunit sequences of the SARS-CoV-2 Spike protein.", style={'fontStyle': 'italic', 'marginBottom': '10px'}),
        html.P("To access the model, please request access on Hugging Face:", style={'marginBottom': '5px'}),
        html.A("https://huggingface.co/SimoRancati/SARITA", href="https://huggingface.co/SimoRancati/SARITA", target="_blank", style={'color': '#1E88E5', 'textDecoration': 'none', 'fontWeight': 'bold'})
    ]),

    # Load model card
    html.Div(style=CARD_STYLE, children=[
        html.Label("Local Model Directory", style=LABEL_STYLE),
        dcc.Input(id="model-dir", type="text", placeholder="/path/to/local/model", style={'width': '100%', 'padding': '8px', 'marginTop': '5px'}),
        html.Button("Load Model", id="load-model", n_clicks=0, style={'marginTop': '15px', 'backgroundColor': '#1E88E5', 'color': '#fff', 'border': 'none', 'padding': '10px 20px', 'borderRadius': '5px', 'cursor': 'pointer'}),
        html.Span(id="status", style={'marginLeft': '15px', 'fontWeight': 'bold'})
    ]),

    html.Hr(),

    # Generation parameters card
    html.Div(style=CARD_STYLE, children=[
        html.Label("Prompt", style=LABEL_STYLE),
        dcc.Input(id="prompt", type="text", placeholder="Starting amino acid sequence (example: MFVFLVLLPLVSSQ, DDFTGCVIAW)", style={'width': '100%', 'padding': '8px', 'marginTop': '5px'}),
        html.Div(style={'display': 'flex', 'justifyContent': 'space-between', 'marginTop': '20px'}, children=[
            html.Div(children=[
                html.Label("Min Length", style=LABEL_STYLE),
                dcc.Input(id="min_len", type="number", value=701, style={'width': '100px', 'padding': '5px'}),
                html.P("The minimum length of the generated sequence.", style={'fontSize': '12px', 'color': '#555'})
            ]),
            html.Div(children=[
                html.Label("Max Length", style=LABEL_STYLE),
                dcc.Input(id="max_len", type="number", value=701, style={'width': '100px', 'padding': '5px'}),
                html.P("The maximum length of the generated sequence.", style={'fontSize': '12px', 'color': '#555'})
            ]),
            html.Div(children=[
                html.Label("Repetition Penalty", style=LABEL_STYLE),
                dcc.Input(id="rep_pen", type="number", step=0.1, value=1.2, style={'width': '100px', 'padding': '5px'}),
                html.P("Penalty for repetition in the generated sequence.", style={'fontSize': '12px', 'color': '#555'})
            ]),
            html.Div(children=[
                html.Label("Sequence Count", style=LABEL_STYLE),
                dcc.Input(id="nseq", type="number", value=1, style={'width': '100px', 'padding': '5px'}),
                html.P("Number of sequences to generate.", style={'fontSize': '12px', 'color': '#555'})
            ])
        ]),

        html.Button("Generate", id="generate", n_clicks=0, disabled=True, style={'marginTop': '20px', 'backgroundColor': '#43A047', 'color': '#fff', 'border': 'none', 'padding': '10px 30px', 'borderRadius': '5px', 'cursor': 'pointer'})
    ]),

    # Output section with spinner
    dcc.Loading(
        id="loading-spinner",
        type="circle",
        children=[
            dcc.Download(id="download-fasta"),
            dcc.Graph(id="aa-plot")
        ],
        style={'textAlign': 'center', 'marginTop': '30px'}
    ),

    # Hidden store
    dcc.Store(id="model-dir-store")
])

@app.callback(
    Output("model-dir-store", "data"),
    Output("status", "children"),
    Output("generate", "disabled"),
    Input("load-model", "n_clicks"),
    State("model-dir", "value"),
    prevent_initial_call=True
)
def load_local_model(n_clicks, model_dir):
    global model, tokenizer
    if not model_dir:
        return None, "❌ Please enter a valid path", True
    try:
        model = AutoModelForCausalLM.from_pretrained(model_dir, trust_remote_code=True)
        tokenizer = AutoTokenizer.from_pretrained(model_dir)
        return model_dir, f"✅ Model loaded from: {model_dir}", False
    except Exception as e:
        return None, f"❌ Load error: {e}", True

@app.callback(
    Output("download-fasta", "data"),
    Output("aa-plot", "figure"),
    Input("generate", "n_clicks"),
    State("prompt", "value"),
    State("min_len", "value"),
    State("max_len", "value"),
    State("rep_pen", "value"),
    State("nseq", "value"),
    prevent_initial_call=True
)
def do_generate(n_clicks, prompt, min_len, max_len, rep_pen, nseq):
    seqs = generate_sequences(
        prompt or "",
        min_len or 701,
        max_len or 701,
        rep_pen or 1.2,
        nseq or 1
    )
    fasta_str = to_fasta(seqs)

    aa_counts = {aa: seqs[0].count(aa) for aa in "ACDEFGHIKLMNPQRSTVWY"}
    fig = px.bar(
        x=list(aa_counts.keys()),
        y=list(aa_counts.values()),
        labels={"x": "Amino Acid", "y": "Count"},
        title="AA Distribution in Sequence 1"
    )
    return dcc.send_string(fasta_str, filename="sarita_sequences.fasta"), fig

if __name__ == "__main__":
    app.run(debug=True)