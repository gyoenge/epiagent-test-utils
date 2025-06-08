import os
import scanpy as sc
import numpy as np
from epiagent.tokenization import tokenization
from epiagent.preprocessing import global_TFIDF
from epiagent.dataset import CellDatasetForDI, collate_fn_for_DI
from torch.utils.data import DataLoader
from epiagent.dataset import CellDataset, collate_fn
from epiagent.model import EpiAgent
import torch
import torch.nn as nn
from epiagent.train import fine_tune_epiagent_for_UFE

def main(
    # finetune parameters
    batch_size = 5, 
    num_steps = 200000,
    learning_rate=1e-4,
    save_steps=20000,
    log_steps=500,
    warmup_steps=10000,
    save_dir = '/root/epiagent/output',
    # input
    data_dir = '/workspace/data',
    model_path = '../model/pretrained_EpiAgent.pth',
): 
    # Save parameters
    os.makedirs(save_dir, exist_ok=True)
    param_txt_path = os.path.join(save_dir, "finetune_params.txt")
    with open(param_txt_path, "w") as f:
        f.write(f"batch_size: {batch_size}\n")
        f.write(f"num_steps: {num_steps}\n")
        f.write(f"learning_rate: {learning_rate}\n")
        f.write(f"save_steps: {save_steps}\n")
        f.write(f"log_steps: {log_steps}\n")
        f.write(f"warmup_steps: {warmup_steps}\n")
        f.write(f"save_dir: {save_dir}\n")
        f.write(f"data_dir: {data_dir}\n")
        f.write(f"model_path: {model_path}\n")

    # Load the raw dataset
    input_path = os.path.join(data_dir, 'sample/raw_h5ad/Kanemaru2023_downsampled_10000_cells.h5ad')
    adata = sc.read_h5ad(input_path)

    # Load the cCRE document frequency data
    cCRE_document_frequency = np.load(os.path.join(data_dir, 'cCRE_document_frequency.npy'))

    # Apply TF-IDF transformation
    global_TFIDF(adata, cCRE_document_frequency)

    # Perform tokenization to create cell sentences
    tokenization(adata)

    # Extract cell sentences from the AnnData object
    cell_sentences = adata.obs['cell_sentences'].tolist()

    # Create the training dataset
    train_cell_dataset = CellDatasetForDI(
        adata=adata,
        cell_sentences=cell_sentences,
        max_length=8192,
        alpha_for_CCA=1,
        num_cCRE=1355445,
        is_random=False
    )

    # Create the training DataLoader
    train_batch_size = batch_size
    train_dataloader = DataLoader(
        train_cell_dataset,
        batch_size=train_batch_size,
        shuffle=True,
        num_workers=16,
        collate_fn=collate_fn_for_DI
    )

    # Set the device (GPU if available)
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

    # Initialize the EpiAgent model with appropriate configurations
    pretrained_model = EpiAgent(
        vocab_size=1355449,
        num_layers=18,
        embedding_dim=512,
        num_attention_heads=8,
        max_rank_embeddings=8192,
        use_flash_attn=True,
        pos_weight_for_RLM=torch.tensor(1.),
        pos_weight_for_CCA=torch.tensor(1.)
    )

    # Load the pre-trained weights into the model
    pretrained_model.load_state_dict(torch.load(model_path))

    # Set criterion for signal reconstruction (SR)
    pretrained_model.criterion_SR = nn.MSELoss()

    # Move the model to the specified device
    pretrained_model.to(device)

    # Fine-tune the model
    fine_tuned_model = fine_tune_epiagent_for_UFE(
        model=pretrained_model,
        train_dataloader=train_dataloader,
        num_steps=num_steps, 
        save_dir=save_dir,
        device=device,
        learning_rate=learning_rate,
        save_steps=save_steps,
        log_steps=log_steps,
        warmup_steps=warmup_steps,
        is_logging=True
    )

    #### test 
    save_path = os.path.join(save_dir, f"finetuned_final.pth")
    torch.save(fine_tuned_model.state_dict(), save_path)

if __name__ == "__main__":
    main(
        # finetune parameters
        batch_size = 5,
        num_steps = 200000,
        learning_rate=1e-4,
        save_steps=20000,
        log_steps=500,
        warmup_steps=10000,
        save_dir = '/workspace/workspace/epiagent/model/finetune_DE_epoch100',
        # input
        data_dir = '/workspace/workspace/epiagent/data',
        model_path = '/workspace/workspace/epiagent/model/pretrained_EpiAgent.pth',
    )


