### image.py ###
# Author: Marc Zepeda
# Date: 2024-12-09

# Import packages
import os
from PIL import Image
import pandas as pd
from ..gen import io

# Image processing methods
def crop(in_dir: str, out_dir: str, box: tuple):
    """
    crop(): crop all photos in a input directory and save them to new directory.
    
    Parameters:
    in_dir (str): Path to the input directory containing images.
    out_dir (str): Path to the destination folder to save cropped images.
    box (tuple): A 4-tuple defining the left, upper, right, and lower pixel coordinates for the crop (left, top, right, bottom).

    Dependencies: os,PIL,io
    """
    # Ensure destination folder exists
    io.mkdir(out_dir)

    # Check box is a tuple with 4 integers
    if all([isinstance(coordinate,int) for coordinate in box])==False:
        KeyError(f"box={box} was not a tuple with 4 integers that define the left, upper, right, and lower pixel coordinates")
    
    # Loop through all files in the source folder
    for filename in os.listdir(in_dir):
        in_path = os.path.join(in_dir, filename)
        out_path = os.path.join(out_dir, filename)
        
        try:
            if filename.lower().endswith(('.png', '.jpg', '.jpeg', '.tiff', '.tif', '.bmp', '.webp')):
                # Open and crop supported image formats
                with Image.open(in_path) as img:
                    cropped_img = img.crop(box)
                    cropped_img.save(out_path)
                    print(f"Successfully processed: {filename}")
        except Exception as e:
            print(f"Error processing {filename}: {e}")

def convert(in_dir: str, out_dir: str, suffix: str):
    """
    convert(): convert image file types of all photos in a input directory and save them to new directory.
    
    Parameters:
    in_dir (str): Path to the folder containing the original images.
    out_dir (str): Path to the folder where converted images will be saved.
    suffix (str): desired image format suffix (e.g., '.png', '.jpeg', '.tiff', '.webp', etc.).

    Dependencies: PIL,os,io
    """
    # Ensure destination folder exists
    io.mkdir(out_dir)
    
    # Loop through all files in the source folder
    for filename in os.listdir(in_dir):
        in_path = os.path.join(in_dir, filename)
        try:
            if filename.lower().endswith(('.png', '.jpg', '.jpeg', '.tiff', '.tif', '.bmp', '.webp')):
                # Open the image
                with Image.open(in_path) as img:
                    # Convert and save the image in the target format
                    base_name, _ = os.path.splitext(filename)
                    destination_path = os.path.join(out_dir, f"{base_name}{suffix.lower()}")
                    img.convert("RGB").save(destination_path, suffix[1:])
                    print(f"Converted {filename} to {suffix[1:]}")
        except Exception as e:
            print(f"Error converting {filename}: {e}")

def combine(in_dir: str, out_dir: str, out_file: str):
    """
    combine(): combine all images in a folder into a single PDF file.
    
    Parameters:
    in_dir (str): Path to the folder containing images.
    out_dir (str): Path for the output PDF file directory.
    out_file (str): PDF filename.

    Dependencies: PIL,os,io
    """
    # Ensure destination folder exists
    io.mkdir(out_dir)

    images = []
    
    # Loop through all files in the folder and collect images
    for filename in sorted(os.listdir(in_dir)):  # Sort files for consistent order
        file_path = os.path.join(in_dir, filename)
        try:
            if filename.lower().endswith(('.png', '.jpg', '.jpeg', '.tiff', '.tif', '.bmp', '.webp')):
                with Image.open(file_path) as img:
                    # Convert image to RGB (required for PDF)
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

# Image information methods
def info(dir: str):
    """
    info(): Extract information from images in a directory as a dataframe.
    
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