#!/bin/bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

streamlit run "$parent_path"/recognition_app.py 
