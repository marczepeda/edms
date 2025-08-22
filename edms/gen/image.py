''' 
Module: image.py
Author: Marc Zepeda
Created: 2024-12-09
Description: Batch image processing

Usage:
[Image processing]
- crop(): crop all photos in a input directory and save them to new directory
- convert(): convert image file types of all photos in a input directory and save them to new directory
- combine(): combine all images in a folder into a single PDF file

[Image information]
- info(): Extract information from images in a directory as a dataframe
'''

# Import packages
import os
from PIL import Image
import pandas as pd
from . import io

# Image processing
def crop(in_dir: str, out_dir: str, box_dims: tuple=None, box_fracs: tuple=None):
    """
    crop(): crop all photos in a input directory and save them to new directory
    
    Parameters:
    in_dir (str): Path to the input directory containing images.
    out_dir (str): Path to the destination folder to save cropped images.
    box_dims (tuple, optional 1): A 4-tuple defining the left, upper, right, and lower absolute pixel coordinates for the crop (left, top, right, bottom).
    box_facs (tuple, optional 1): A 4-tuple defining the left, upper, right, and lower fraction pixel coordinates for the crop (left, top, right, bottom).

    Dependencies: os,PIL,io
    """
    # Ensure output folder exists
    io.mkdir(out_dir)

    if box_dims is not None:
        # Check box_dims is a tuple with 4 integers 
        if all([isinstance(coordinate,int) for coordinate in box_dims])==False and len(box_dims)==4:
            KeyError(f"box_dims={box_dims} was not a tuple with 4 integers that define the left, upper, right, and lower absolute pixel coordinates")

        # Loop through all files in the input folder
        for filename in os.listdir(in_dir):
            in_path = os.path.join(in_dir, filename)
            out_path = os.path.join(out_dir, filename)
            
            try: # Open and crop supported image formats
                if filename.lower().endswith(('.png', '.jpg', '.jpeg', '.tiff', '.tif', '.bmp', '.webp')):
                    with Image.open(in_path) as img:
                        cropped_img = img.crop(box_dims)
                        cropped_img.save(out_path)
                        print(f"Successfully processed: {filename}")
            except Exception as e:
                print(f"Error processing {filename}: {e}")
    
    elif box_fracs is not None:
        # Check box_fracs is a tuple with 4 floats
        if all([isinstance(coordinate,float) for coordinate in box_fracs])==False and len(box_fracs)==4:
            KeyError(f"box_fracs={box_fracs} was not a tuple with 4 floats that define the left, upper, right, and lower fraction pixel coordinates")

        # Loop through all files in the input folder
        for filename in os.listdir(in_dir):
            in_path = os.path.join(in_dir, filename)
            out_path = os.path.join(out_dir, filename)
            
            try: # Open and crop supported image formats
                if filename.lower().endswith(('.png', '.jpg', '.jpeg', '.tiff', '.tif', '.bmp', '.webp')):
                    with Image.open(in_path) as img:
                        box_dims = (box_fracs[0]*img.width,box_fracs[1]*img.height,box_fracs[2]*img.width,box_fracs[3]*img.height)
                        cropped_img = img.crop(box_dims)
                        cropped_img.save(out_path)
                        print(f"Successfully processed: {filename}")
            except Exception as e:
                print(f"Error processing {filename}: {e}")
         
    else:
        KeyError(f"Neither box_dims nor box_fracs was provided")

def convert(in_dir: str, out_dir: str, suffix: str):
    """
    convert(): convert image file types of all photos in a input directory and save them to new directory
    
    Parameters:
    in_dir (str): Path to the folder containing the original images.
    out_dir (str): Path to the folder where converted images will be saved.
    suffix (str): desired image format suffix (e.g., '.png', '.jpeg', '.tiff', '.webp', etc.).

    Dependencies: PIL,os,io
    """
    # Ensure output folder exists
    io.mkdir(out_dir)
    
    # Loop through all files in the input folder
    for filename in os.listdir(in_dir):
        in_path = os.path.join(in_dir, filename)

        try: # Open the image
            if filename.lower().endswith(('.png', '.jpg', '.jpeg', '.tiff', '.tif', '.bmp', '.webp')):
                with Image.open(in_path) as img: # Convert and save the image in the target format
                    base_name, _ = os.path.splitext(filename)
                    out_path = os.path.join(out_dir, f"{base_name}{suffix.lower()}")
                    img.convert("RGB").save(out_path, suffix[1:])
                    print(f"Converted {filename} to {suffix[1:]}")
        except Exception as e:
            print(f"Error converting {filename}: {e}")

def combine(in_dir: str, out_dir: str, out_file: str):
    """
    combine(): combine all images in a folder into a single PDF file
    
    Parameters:
    in_dir (str): Path to the folder containing images.
    out_dir (str): Path for the output PDF file directory.
    out_file (str): PDF filename.

    Dependencies: PIL,os,io
    """
    # Ensure output folder exists
    io.mkdir(out_dir)

    images = []
    
    # Loop through all files in the folder and collect images
    for filename in sorted(os.listdir(in_dir)):  # Sort files for consistent order
        file_path = os.path.join(in_dir, filename)
        try:
            if filename.lower().endswith(('.png', '.jpg', '.jpeg', '.tiff', '.tif', '.bmp', '.webp')):
                with Image.open(file_path) as img: # Convert image to RGB (required for PDF)
                    images.append(img.convert("RGB"))
        except Exception as e:
            print(f"Skipping {filename}: {e}")
    
    if not images:
        print("No valid images found to combine.")
        return
    
    # Save the images as a single PDF
    try:
        # The first image is used as the primary document; others are appended
        images[0].save(os.path.join(out_dir,out_file), save_all=True, append_images=images[1:])
        print(f"{out_file} saved in {out_dir}")
    except Exception as e:
        print(f"Error creating PDF: {e}")

# Image information
def info(dir: str) -> pd.DataFrame:
    """
    info(): Extract information from images in a directory as a dataframe
    
    Parameters:
    dir (str): Path to the directory containing images.
        
    Depedencies: PIL,os,pandas
    """
    image_info_list = []
    
    for filename in os.listdir(dir):
        file_path = os.path.join(dir, filename)
        try:
            with Image.open(file_path) as img:
                # Extract image details
                image_info = {
                    "Filename": filename,
                    "Width": img.width,
                    "Height": img.height,
                    "Mode": img.mode,
                    "Format": img.format,
                    "Size (bytes)": os.path.getsize(file_path)
                }
                image_info_list.append(image_info)
        except Exception as e:
            print(f"Skipping {filename}: {e}")
    
    return pd.DataFrame(image_info_list)