#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import streamlit as st
from PIL import Image 
import numpy as np 
import glob 
import streamlit.components.v1 as components
import os 
import json 
import argparse
from pathlib import Path


def main(args):
    """
        A simple application to easily classify images with keyboard shortcut.

        The dataset to annotate respect the following structure:
            args.base_path/images/*.png
            args.base_path/classes/*.json

        Pre-annotations can be imported by placing files in args.base_path/classes/.
        Annotation files respect the following schema:
        {
            "label": The image label (Example: "markush_structure")
        }
    """
    st.set_page_config(layout="wide")

    if "skip_annotate" not in st.session_state:
        st.session_state.skip_annotated = False
    
    if "index" not in st.session_state:
        st.session_state.index = "0"

    images_paths = [p for p in glob.glob(args.base_path + "/images/**", recursive=True) if ".png" in p]
    if images_paths == []:
        print("The folder containing images to annotate is empty.")
        exit(0)
        
    # Print number of images annotated and the first image which is not annotated.
    first_image_printed = False
    nb_annotated_images = 0
    for p in glob.glob(args.base_path + "/images/**", recursive=True):
        if not(".png" in p):
            continue
        parts = Path(p).parts
        image_filename = str(Path(*parts[parts.index(("/images").split("/")[-1]) + 1:]).stem)
        
        if not(os.path.exists(args.base_path + f"/classes/{image_filename}.json")):
            if not(first_image_printed):
                print(f"First image which is not annotated: {p}")
                first_image_printed = True
            continue
        nb_annotated_images += 1
    print(f"Number of images annotated: {nb_annotated_images}/{len(images_paths)}")

    # Set callbacks
    def increment_callback():
        st.session_state.index = str(int(st.session_state.index) + 1)
        update_filename()

    def decrement_callback():
        st.session_state.index = str(int(st.session_state.index) - 1)
        update_filename()

    def standard_molecule_callback():
        with open(st.session_state.manual_annotation_path, "w") as f:
            d = {"label": "standard_molecule"}
            json.dump(d, f)
        st.session_state.index = str(int(st.session_state.index) + 1)
        update_filename()

    def markush_structure_callback():
        with open(st.session_state.manual_annotation_path, "w") as f:
            d = {"label": "markush_structure"}
            json.dump(d, f)
        st.session_state.index = str(int(st.session_state.index) + 1)
        update_filename()

    def miscellaneous_callback():
        with open(st.session_state.manual_annotation_path, "w") as f:
            d = {"label": "miscellaneous"}
            json.dump(d, f)
        st.session_state.index = str(int(st.session_state.index) + 1)
        update_filename()

    def segmenation_annotation_error_callback():
        with open(st.session_state.manual_annotation_path, "w") as f:
            d = {"label": "segmentation_annotation_error"}
            json.dump(d, f)
        st.session_state.index = str(int(st.session_state.index) + 1)
        update_filename()

    def polymer_callback():
        with open(st.session_state.manual_annotation_path, "w") as f:
            d = {"label": "polymer"}
            json.dump(d, f)
        st.session_state.index = str(int(st.session_state.index) + 1)
        update_filename()

    def update_filename():
        st.session_state.current_filename = images_paths[int(st.session_state.index)]

    def update_index():
        st.session_state.index = str(images_paths.index(st.session_state.current_filename))

    if "current_filename" not in st.session_state:
        update_filename()
    
    # Set buttons
    zero_col, one_col, two_col, three_col, four_col, five_col, six_col, _ = st.columns([1, 1, 1, 1, 2, 1, 1, 1])
    with zero_col:
        st.button('increment_(0)' , on_click=increment_callback)
    with one_col:
        st.button('standard_molecule_(1)' , on_click=standard_molecule_callback)
    with two_col:
        st.button('markush_structure_(2)' , on_click=markush_structure_callback)
    with three_col:
        st.button('miscellaneous_(3)' , on_click=miscellaneous_callback)
    with four_col:
        st.button('segmentation_annotation_error_(4)' , on_click=segmenation_annotation_error_callback)
    with five_col:
        st.button('polymer_(5)' , on_click=polymer_callback)
    with six_col:
        st.button('decrement_(6)' , on_click=decrement_callback)    

    parts = Path(images_paths[int(st.session_state.index)]).parts
    image_filename = str(Path(*parts[parts.index(("/images").split("/")[-1]) + 1:]).stem)
    
    # Display Image Filename 
    st.text_input("Dataset Filename:", key="current_filename", on_change=update_index)

    # Display Manual Annotation
    st.session_state.manual_annotation_path = args.base_path + f"/classes/{image_filename}.json"
    if os.path.exists(st.session_state.manual_annotation_path):
        with open(st.session_state.manual_annotation_path) as f:
            manual_annotation = json.load(f)
            st.markdown(f"Manual annotation: {manual_annotation['label']}", unsafe_allow_html=True)  

    # Display Image           
    if st.session_state.index != "":
        image_path = images_paths[int(st.session_state.index)]
        image = Image.open(image_path)
        image = np.array(image)
        st.image(image)

    # Use Keyboard Shortcuts
    components.html(
        """
        <script>
        const doc = window.parent.document;
        buttons = Array.from(doc.querySelectorAll('button'));
        const zero_button = buttons.find(el => el.innerText === 'increment_(0)');
        const one_button = buttons.find(el => el.innerText === 'standard_molecule_(1)');
        const two_button = buttons.find(el => el.innerText === 'markush_structure_(2)');
        const three_button = buttons.find(el => el.innerText === 'miscellaneous_(3)');
        const four_button = buttons.find(el => el.innerText === 'segmentation_annotation_error_(4)');
        const five_button = buttons.find(el => el.innerText === 'polymer_(5)');
        const six_button = buttons.find(el => el.innerText === 'decrement_(6)');
        doc.addEventListener('keydown', function(e) {
            switch (e.keyCode) {
                case 96: 
                    zero_button.click();
                    break;
                case 97: 
                    one_button.click();
                    break;
                case 98: 
                    two_button.click();
                    break;
                case 99:
                    three_button.click();
                    break;
                case 100:
                    four_button.click();
                    break;
                case 101:
                    five_button.click();
                    break;
                case 102:
                    six_button.click();
                    break;
            }
        });
        </script>
        """,
        height=0,
        width=0,
    )
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--base-path', type = str, default="../../data/classification/default")
    args = parser.parse_args()
    main(args)

