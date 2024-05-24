# MolAnnotator

This repository implements applications to ease manual annotation of chemical images. 

## Installation

```
conda create -n molannotator python=3.11
conda activate molannotator
pip install -e .
```

## Classification
The classification application allows to assign to an image one of the following categories: `standard_molecule` (Molecular Structure), `markush_structure` (Markush Structure), `miscellaneous` (Background), `polymer` (Markush Structure), and `segmentation_annotation_error`.

To run the application:
1. Move the images to be annotated to `./data/classfication/default/images/*.png`,
2. (optional) Move pre-annotations to `./data/classfication/default/classes/*.json`,
3. Run `bash ./molannotator/apps/run_classification_app.sh`.

Annotations are saved in `./data/classfication/default/classes/*.json`.

## Recognition
The recognition application allows to annotate the molecular graph of a chemical-structure image.

To run the application:
1. Move the images to be annotated to `./data/recognition/default/images/*.png`,
2. (optional) Create pre-annotations by running `python3 ./molannotator/molscribe_pre_annotate.py` (see [MolScribe](https://github.com/thomas0809/MolScribe) for installation steps),
3. (optional) Move pre-annotations to `./data/recognition/default/predictions/*.mol`,
3. Run `bash ./molannotator/apps/run_recognition_app.sh`.

Annotations are saved in `./data/recognition/default/molfiles/*.mol`.
