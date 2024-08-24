# This script will read the images from the designated folder, apply blur + grey scale + contour and calculate the average cell size and average cell area
# the threshold can be treated according to different cell size/ roundness/ area
# output image: original image + processed image + extracted cell
import cv2
import numpy as np
import matplotlib.pyplot as plt
import os

image_folder = 'F:/MScR Reg Med + Tissue Repair/02 - Project Two/00-Analysis/Videos/pytest/Laminar-1/'

def process_image(image_path):
    image = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)
    blurred = cv2.GaussianBlur(image, (11, 11), 0)
    thresh = cv2.adaptiveThreshold(blurred, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY_INV, 11, 2)
    contours, _ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    
    cell_areas = []
    mask = np.zeros_like(image)
    for contour in contours:
        if len(contour) < 5:  # fitEllipse requires at least 5 points
            continue
        
        ellipse = cv2.fitEllipse(contour)
        (x, y), (MA, ma), angle = ellipse
        a = MA / 2
        b = ma / 2
        area = np.pi * a * b
        perimeter = cv2.arcLength(contour, True)
        
        if perimeter == 0:
            continue
        
        circularity = 4 * np.pi * (area / (perimeter * perimeter))
        
        # Debug print statements
        print(f'Contour properties - Area: {area:.2f}, Circularity: {circularity:.2f}, Semi-major axis: {a:.2f}')
        
        # Check if the cell satisfies the conditions
        if 200 < area < 1500 and 0.3 < circularity < 1.5 and 14 < a < 30:
            cv2.ellipse(mask, ellipse, (255), thickness=-1)
            cell_areas.append(area)
    
    extracted_cells = cv2.bitwise_and(image, image, mask=mask)
    return len(cell_areas), cell_areas, image, thresh, extracted_cells

def save_images(image_list, output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    for idx, (filename, _, _, original, processed, extracted) in enumerate(image_list):
        base_filename = os.path.splitext(filename)[0]
        cv2.imwrite(os.path.join(output_folder, f'{base_filename}_original.png'), original)
        cv2.imwrite(os.path.join(output_folder, f'{base_filename}_processed.png'), processed)
        cv2.imwrite(os.path.join(output_folder, f'{base_filename}_extracted.png'), extracted)


output_folder = os.path.join(image_folder, 'output_images')

results = []

for filename in os.listdir(image_folder):
    if filename.lower().endswith(('.png', '.jpg', '.jpeg', '.tiff', '.bmp', '.gif')):
        image_path = os.path.join(image_folder, filename)
        num_cells, cell_areas, original, processed, extracted = process_image(image_path)
        results.append((filename, num_cells, cell_areas, original, processed, extracted))

save_images(results, output_folder)

total_cells = [num_cells for _, num_cells, _, _, _, _ in results]
all_cell_areas = [area for _, _, cell_areas, _, _, _ in results for area in cell_areas]

average_num_cells = np.mean(total_cells)
std_num_cells = np.std(total_cells)

average_area = np.mean(all_cell_areas)
std_area = np.std(all_cell_areas)

print(f'Average Number of Cells: {average_num_cells} ± {std_num_cells}')
print(f'Average Area of Cells: {average_area:.2f} ± {std_area:.2f} pixel^2')

if results:
    plt.figure(figsize=(18, 6))

    plt.subplot(1, 3, 1)
    plt.title('Original Image')
    plt.imshow(results[0][3], cmap='gray')

    plt.subplot(1, 3, 2)
    plt.title('Processed Image')
    plt.imshow(results[0][4], cmap='gray')

    plt.subplot(1, 3, 3)
    plt.title('Extracted Cells')
    plt.imshow(results[0][5], cmap='gray')

    plt.tight_layout()
    plt.show()

plt.figure(figsize=(12, 6))

plt.subplot(1, 2, 1)
plt.title('Average Cell Areas')
plt.hist(all_cell_areas, bins=20, color='blue', edgecolor='black')
plt.axvline(average_area, color='red', linestyle='dashed', linewidth=1)
plt.axvline(average_area - std_area, color='green', linestyle='dashed', linewidth=1)
plt.axvline(average_area + std_area, color='green', linestyle='dashed', linewidth=1)
plt.xlabel('Area (pixel^2)')
plt.ylabel('Frequency')

plt.subplot(1, 2, 2)
plt.title('Average Number of Cells')
plt.bar(['Average'], [average_num_cells], yerr=[std_num_cells], color='green', capsize=10)
plt.ylabel('Number of Cells')

plt.tight_layout()
plt.show()