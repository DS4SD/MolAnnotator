#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import torch
from molscribe import MolScribe
from huggingface_hub import hf_hub_download
import glob 
from rdkit import Chem
from rdkit.Chem import rdmolfiles
from tqdm import tqdm 
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--base-path', type = str, default="../../data/recognition/default")
    args = parser.parse_args()

    images_paths = [p for p in glob.glob(f"{args.base_path}/rec/images/**/*.png", recursive=True)]
    
    ckpt_path = hf_hub_download('yujieq/MolScribe', 'swin_base_char_aux_1m.pth')
    model = MolScribe(ckpt_path, device=torch.device('cuda'))
    
    for image_path in tqdm(images_paths):
        output = model.predict_image_file(image_path, return_atoms_bonds=True, return_confidence=True)
        molfile_block = output["molfile"]
        molecule = rdmolfiles.MolFromMolBlock(molfile_block)
        if molecule is None:
            molecule = Chem.MolFromSmiles("C")
            
        molfile_path = f"{args.base_path}/rec/predictions/" + "/".join(image_path.split("/")[-2:])[:-4] + ".mol"
        rdmolfiles.MolToMolFile(
            molecule, 
            molfile_path, 
            kekulize = False 
        )
        
if __name__ == "__main__":
    main()
    