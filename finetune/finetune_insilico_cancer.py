import os
import scanpy as sc
from epiagent.dataset import CellDatasetForUFE, collate_fn_for_UFE
from torch.utils.data import DataLoader
from epiagent.model import EpiAgent
import torch
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
    data_path = '/workspace/epiagent/data/insilico/processed_h5ad/downsampled_train_20000.h5ad',
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
        f.write(f"data_path: {data_path}\n")
        f.write(f"model_path: {model_path}\n")

    # Load the preprocessed AnnData
    input_path = os.path.join(data_path)
    adata = sc.read_h5ad(input_path)

    # Extract cell sentences from the AnnData object
    cell_sentences = adata.obs['cell_sentences'].tolist()

    # Create the training dataset
    train_cell_dataset = CellDatasetForUFE(
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
        collate_fn=collate_fn_for_UFE
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

    # Ensure the CCA loss uses a positive weight of 1
    pretrained_model.criterion_CCA.pos_weight = torch.tensor(1.)

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
    # main(
    #     # finetune parameters
    #     batch_size = 5,
    #     num_steps = 200000,
    #     learning_rate=1e-4,
    #     save_steps=20000,
    #     log_steps=500,
    #     warmup_steps=10000,
    #     save_dir = '/workspace/epiagent/model/finetune_IS_epoch50',
    #     # input
    #     data_path = '/workspace/epiagent/data/insilico/processed_h5ad/downsampled_train_20000_cellsentenced.h5ad',
    #     model_path = '/workspace/epiagent/model/pretrained_EpiAgent.pth',
    # )
    main(
        # finetune parameters
        batch_size = 5,  # 1 epoch = 3571 step 
        num_steps = 178550,  # 50 epochs = 178550 steps
        learning_rate=1e-4,
        save_steps=17855,  # save every 5 epochs
        log_steps=714,  # log every 5 times per 1 epoch
        warmup_steps=10000,  
        save_dir = '/workspace/epiagent/model/finetune_IS_epoch50_stable_17855',
        # input
        data_path = '/workspace/epiagent/data/insilico/processed_h5ad/downsampled_stable_17855_cellsentenced.h5ad',
        model_path = '/workspace/epiagent/model/pretrained_EpiAgent.pth',
    )