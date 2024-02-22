#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import streamlit as st
from rdkit import Chem 
from streamlit_ketcher import st_ketcher
from PIL import Image 
import numpy as np 
import glob 
import os 
import random 
import argparse
from pathlib import Path

from mol_depict.utils.utils_drawing import draw_molecule_rdkit


def main(args):
    """
        A simple application to easily annotate molecular graphs from chemical-structure images.

        The dataset to annotate respect the following structure:
            args.base_path/images/*.png
            args.base_path/molfiles/*.mol
            args.base_path/predictions/*.mol

        Pre-annotations can be imported by placing files in args.base_path/predictions/, and obtained using the script molannotator/molscribe_pre_annotate.py.
    """
    # Set layout
    st.set_page_config(layout="wide")

    # Set session state attributes
    if "index" not in st.session_state:
        st.session_state.index = "0"

    if "molecule" not in st.session_state:
        st.session_state.molecule = None

    if "sub_index" not in st.session_state:
        st.session_state.sub_index = "0"
    
    # Set images, predictions and annotations paths
    predictions_molfiles_paths = sorted([p for p in glob.glob(args.base_path + "/predictions/**.mol", recursive=True)])
    images_paths = sorted([p for p in glob.glob(args.base_path + "/images/**.png", recursive=True)])
    annotation_molfiles_paths = sorted([args.base_path + "/molfiles/" + p.split("/")[-1][:-4] + ".mol" for p in images_paths])

    print("Number of predictions: ", len(predictions_molfiles_paths))
    print("Number of images: ", len(images_paths))
    
    # Set callbacks
    def increment_sub_index():
        st.session_state.sub_index = str(int(st.session_state.sub_index) + 1)
        st.session_state.index = str(subindices[int(st.session_state.sub_index)])
        update_filename()

    def decrement_sub_index():
        st.session_state.sub_index = str(int(st.session_state.sub_index) - 1)
        st.session_state.index = str(subindices[int(st.session_state.sub_index)])
        update_filename()

    def select_random_index():
        to_annotate_molfiles_index = [i for i, p in enumerate(annotation_molfiles_paths) if not os.path.exists(p)]
        st.session_state.index = str(random.choice(to_annotate_molfiles_index))
        update_filename()

    def increment_index():
        st.session_state.index = str(int(st.session_state.index) + 1)
        update_filename()
    
    def decrement_index():
        st.session_state.index = str(int(st.session_state.index) - 1)
        update_filename()

    def save_molecule():
        Chem.MolToMolFile(
            st.session_state.molecule,
            annotation_molfiles_paths[int(st.session_state.index)]
        )

    def update_filename():
        st.session_state.current_filename = images_paths[int(st.session_state.index)].split("/")[-1]

    def update_index():
        st.session_state.index = str(images_paths.index(args.base_path + "/images/" + st.session_state.current_filename))
        update_filename()

    if "current_filename" not in st.session_state:
        update_filename()

    # Set buttons
    col1, col2, col3, col4, col5, col6, _ = st.columns([1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8])
    with col1:
        st.button('Save annotation', on_click = save_molecule)
    with col2:
        st.button('Select random', on_click = select_random_index)
    with col3:
        st.button('Increment', on_click = increment_index)
    with col4:
        st.button('Decrement', on_click = decrement_index)
    with col5:
        st.button('Increment (not-annotated files)', on_click = increment_sub_index)
    with col6:
        st.button('Decrement (not-annotated files)', on_click = decrement_sub_index)
    
    if os.path.exists(annotation_molfiles_paths[int(st.session_state.index)]):
        st.markdown(f"<font color='green'> Annotated </font>", unsafe_allow_html=True) 
    else:
        st.markdown(f"<font color='red'> Not Annotated </font>", unsafe_allow_html=True) 

    st.text_input("Dataset Index:", key="index", on_change=update_filename)
    st.text_input("Dataset Filename:", key="current_filename", on_change=update_index)
    col1, col2 = st.columns(2)

    if st.session_state.index != "":
        # Read Image
        image_path = images_paths[int(st.session_state.index)]
        image = Image.open(image_path)
        image = np.array(image)

        # Read Molecule
        if os.path.exists(args.base_path + "/molfiles/" + st.session_state.current_filename[:-4] + ".mol"):
            molfile_path = args.base_path + "/molfiles/" + st.session_state.current_filename[:-4] + ".mol"
            molecule = Chem.MolFromMolFile(molfile_path, sanitize=False, removeHs=False, strictParsing=False)
        elif os.path.exists(args.base_path + "/predictions/" + st.session_state.current_filename[:-4] + ".mol"):
            molfile_path = args.base_path + "/predictions/" + st.session_state.current_filename[:-4] + ".mol"
            molecule = Chem.MolFromMolFile(molfile_path, sanitize=False, removeHs=False, strictParsing=False)
        else:
            molfile_path = ""
            molecule = Chem.MolFromSmiles("")
        input_molblock = Chem.rdmolfiles.MolToMolBlock(molecule, kekulize=False)
        
        with col1:
            annotated_molfile = st_ketcher(input_molblock, molecule_format="MOLFILE")
            st.session_state.molecule = Chem.MolFromMolBlock(annotated_molfile, sanitize=False, removeHs=False, strictParsing=False)
            
            if st.session_state.molecule:
                # Display
                smiles = Chem.MolToSmiles(st.session_state.molecule)
                st.markdown(f"Image path: {image_path}")
                st.markdown(f"Read molfile path: {molfile_path}")
                st.markdown(f"Save molfile path: {args.base_path}/molfiles/{st.session_state.current_filename[:-4]}.mol")
                st.markdown(f"SMILES: ``{smiles}``")
                annotated_molfile_display = annotated_molfile.replace('\n', '<br>')
                st.markdown(f"Molfile: {annotated_molfile_display}", unsafe_allow_html=True)

        with col2:
            st.image(image)

    st.markdown(f"RDKit Image: ")
    rdkit_image = draw_molecule_rdkit(molecule = molecule, smiles = smiles, augmentations = False, display_atom_indices = True)
    rdkit_image = np.array(rdkit_image.permute(1, 2, 0))/255
    st.image(rdkit_image)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--base-path', type = str, default="../../data/recognition/default")
    args = parser.parse_args()

    images_paths = sorted([p for p in glob.glob(args.base_path + "/images/**.png", recursive=True)])
    if images_paths == []:
        print("The folder containing images to annotate is empty.")
        exit(0)
        
    # Get not annotated indices
    s = 0
    subindices = []
    for image_path in images_paths:
        parts = Path(image_path).parts
        filename = str(Path(*parts[parts.index(("/images").split("/")[-1]) + 1:]).stem)

        if not os.path.exists(args.base_path + "/molfiles/" + filename + ".mol"):
            s += 1
            subindices.append(images_paths.index(args.base_path + "/images/" + filename + ".png"))
            
    print("Number of annotations: ", len(images_paths) - s)
    main(args)