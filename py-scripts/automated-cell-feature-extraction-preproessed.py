import cv2
import numpy as np
import matplotlib.pyplot as plt
import os

# Define a function to process a single image and return cell information
def process_image(image_path):
    image = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)
    blurred = cv2.GaussianBlur(image, (11, 11), 0)
    thresh = cv2.adaptiveThreshold(blurred, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY_INV, 11, 2)
    contours, _ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    
    cell_areas = []
    mask = np.zeros_like(image)
    for contour in contours:
        area = cv2.contourArea(contour)
        perimeter = cv2.arcLength(contour, True)
        if perimeter == 0:
            continue
        circularity = 4 * np.pi * (area / (perimeter * perimeter))
        (x, y), radius = cv2.minEnclosingCircle(contour)
        if 50 < area < 1500 and 0.7 < circularity < 1.2 and 5 < radius < 30:
            cv2.circle(mask, (int(x), int(y)), int(radius), (255), thickness=-1)
            cell_areas.append(area)
    
    extracted_cells = cv2.bitwise_and(image, image, mask=mask)
    return len(cell_areas), cell_areas, extracted_cells

# Simulating the process for 10 images using the provided image
image_path = 'F:/MScR Reg Med + Tissue Repair/02 - Project Two/00-Analysis/Videos/2024-06-25_00-21-54-515-edit-substract-background.jpg'
results = []
for i in range(10):
    num_cells, cell_areas, extracted_cells = process_image(image_path)
    results.append((f'Image_{i+1}', num_cells, cell_areas))

# Print the results
for image_file, num_cells, cell_areas in results:
    print(f'Image: {image_file}')
    print(f'Number of Cells: {num_cells}')
    print(f'Areas of Cells: {cell_areas}')
    print('---')
