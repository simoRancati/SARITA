# SARITA
SARITA, a GPT-3-based LLM with 1.2B parameters, generates synthetic SARS-CoV-2 Spike S1 sequences. Built on the RITA protein model, it uses continuous learning. When trained on Alpha, Beta, and Gamma variants (data up to Feb 2021), SARITA accurately predicts future S1 mutations

This model is a continuous learning of [lightonai/RITA_s](https://huggingface.co/lightonai/RITA_s).
It achieves the following results on the evaluation set:
- Loss: 0.0293
- Accuracy: 0.9914


Model | #Params | d_model | layers 
--- | --- | --- | --- | 
[Small](https://huggingface.co/SimoRancati/SARITA_S) | 85M  | 768 | 12 
[Medium](https://huggingface.co/SimoRancati/SARITA_M) | 300M | 1024 | 24 
[Large](https://huggingface.co/SimoRancati/SARITA_L)| 680M | 1536 | 24 
[XLarge](https://huggingface.co/SimoRancati/SARITA_XL)| 1.2B | 2048 | 24 


## Model description

SARITA S is an LLM with 85 million parameters, based on GPT-3 architecture, designed to generate high-quality synthetic SARS-CoV-2 Spike sequences. 
SARITA is trained via continuous learning on the pre-existing protein model RITA.

## Intended uses & limitations

This model can be used by user to generate synthetic Spike proteins of SARS-CoV-2 Virus. 


### Training hyperparameters

The following hyperparameters were used during training:
- learning_rate: 5e-05
- train_batch_size: 48
- eval_batch_size: 48
- seed: 42
- optimizer: Adam with betas=(0.9,0.999) and epsilon=1e-08
- lr_scheduler_type: linear
- num_epochs: 10
- mixed_precision_training: Native AMP

### Training results

| Training Loss | Epoch | Step   | Validation Loss | Accuracy |
|:-------------:|:-----:|:------:|:---------------:|:--------:|
| 0.0303        | 1.0   | 10013  | 0.0302          | 0.9912   |
| 0.0297        | 2.0   | 20026  | 0.0300          | 0.9912   |
| 0.0294        | 3.0   | 30039  | 0.0295          | 0.9913   |
| 0.0294        | 4.0   | 40052  | 0.0294          | 0.9913   |
| 0.0293        | 5.0   | 50065  | 0.0294          | 0.9913   |
| 0.0292        | 6.0   | 60078  | 0.0293          | 0.9914   |
| 0.029         | 7.0   | 70091  | 0.0293          | 0.9914   |
| 0.0288        | 8.0   | 80104  | 0.0293          | 0.9914   |
| 0.0285        | 9.0   | 90117  | 0.0294          | 0.9914   |
| 0.0283        | 10.0  | 100130 | 0.0295          | 0.9914   |


### Framework versions

- Transformers 4.20.1
- Pytorch 1.9.0+cu111
- Datasets 2.18.0
- Tokenizers 0.12.1


### Usage 
Instantiate a model like so:
``` python
from transformers import AutoModelForCausalLM, AutoTokenizer
model = AutoModelForCausalLM.from_pretrained("SimoRancati/SARITA_*", trust_remote_code=True)
tokenizer = AutoTokenizer.from_pretrained("SimoRancati/SARITA_*")
```
for generation used this code:
``` python
# Check for GPU availability and move the model to GPU
device = "cuda" if torch.cuda.is_available() else "cpu"
model = model.to(device)

start = ['MFVFLVLLPLVSSQ']

for i in range(len(start)):
    # Prepare model inputs
    model_inputs = tokenizer([start[i]], return_tensors="pt")
    model_inputs = {k: v.to(device) for k, v in model_inputs.items()}

    # Generate predictions using the model
    generated_ids = model.generate(**model_inputs, min_length=701, max_length=701,
                                   do_sample=True, top_k=950, repetition_penalty=1.2,
                                   num_return_sequences=100, eos_token_id=2, truncation=True)

    # Decode and print outputs
    generated_sequences = []
    for f in range(len(generated_ids)):
        sequence = tokenizer.decode(generated_ids[f], skip_special_tokens=True).replace(' ', '')
        generated_sequences.append(sequence)
```
# Repository
This Repository contains all the code used to built the [datatset](/SimoRancati/SARITA/Dataset), to [train the model](/SimoRancati/SARITA/Training), to [generate synthetic sequences](/SimoRancati/SARITA/Generation) and to [evaluate the different models](/SimoRancati/SARITA/Evaluation)

# IMPORTANT 
SARITA model is public, but downloading it requires approval.  
To request access, go to [HuggingFace]( https://huggingface.co/SimoRancati) and click on the **Request Access** button and provide a brief explanation of your intended use.

## License
The use of this model is restricted to research purposes. Commercial use is not allowed without prior approval.
