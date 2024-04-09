import os
import cv2
import numpy as np
from skimage import filters, color, io, exposure

# Function to calculate contrast (standard deviation of pixel values)
def calculate_contrast(image):
    return np.std(image)

# Function to calculate luminance
def calculate_luminance(image):
    # Convert to YUV and take the Y channel as an approximation of luminance
    yuv_image = cv2.cvtColor(image, cv2.COLOR_BGR2YUV)
    y_channel = yuv_image[:, :, 0]
    return np.average(y_channel)

# Function to calculate visual complexity (edge density)
def calculate_complexity(image):
    gray_image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    edges = filters.sobel(gray_image)
    edge_density = np.sum(edges > np.mean(edges)) / np.prod(edges.shape)
    return edge_density

# Function to adjust the luminance of an image
def adjust_luminance(image, target_luminance, white_threshold=245):
    # Convert to YUV and take the Y channel as an approximation of luminance
    yuv_image = cv2.cvtColor(image, cv2.COLOR_BGR2YUV)
    y_channel, u_channel, v_channel = cv2.split(yuv_image)

    # Determine which pixels are not white
    not_white_mask = np.all(image < white_threshold, axis=2)

    # Calculate the average luminance of non-white pixels
    current_luminance = np.mean(y_channel[not_white_mask])

    # Calculate the current luminance and the necessary scaling factor
    current_luminance = np.mean(y_channel)
    scaling_factor = target_luminance / current_luminance

    # Scale the Y channel
    y_channel = np.clip(y_channel * scaling_factor, 0, 255).astype(y_channel.dtype)

    # Merge the channels back and convert to BGR
    adjusted_yuv = cv2.merge((y_channel, u_channel, v_channel))
    adjusted_img = cv2.cvtColor(adjusted_yuv, cv2.COLOR_YUV2BGR)

    return adjusted_img

# Function to adjust the contrast of an image
def adjust_contrast(image, target_std):
    img_std = np.std(image)
    if img_std > 0:
        adjusted_img = exposure.rescale_intensity(image, in_range=(np.percentile(image, 2), np.percentile(image, 98)))
        adjusted_img = ((adjusted_img - np.mean(adjusted_img)) / np.std(adjusted_img)) * target_std + np.mean(adjusted_img)
    else:
        adjusted_img = image
    return np.clip(adjusted_img, 0, 255).astype(np.uint8)

# Load images
dir_path = os.getcwd()
file_path = os.path.join(dir_path, "stimuli", "originals")
files = [os.path.join("stimuli", "originals", item) for item in os.listdir(file_path) if item.endswith(".png")]
file_names = [os.path.join("stimuli", item) for item in os.listdir(file_path) if item.endswith(".png")]
files.sort()
file_names.sort()

images = [cv2.imread(file) for file in files]

# Calculate metrics for each image
contrasts = [calculate_contrast(image) for image in images]
luminances = [calculate_luminance(image) for image in images]
complexities = [calculate_complexity(image) for image in images]

mean_contrast = np.mean(contrasts)
std_contrast = np.std(contrasts)
mean_luminance = np.mean(luminances)
std_luminance = np.std(luminances)

# Normalize images
normalized_images = []
for image in images:
    # image = images[0]

    # Adjust luminance
    image_luminance = calculate_luminance(image)
    if image_luminance < mean_luminance - std_luminance:
        target_luminance = mean_luminance - std_luminance
        image = adjust_luminance(image, target_luminance)
    elif image_luminance > mean_luminance + std_luminance:
        target_luminance = mean_luminance + std_luminance
        image = adjust_luminance(image, target_luminance)

    # Adjust contrast
    image_contrast = calculate_contrast(image)
    if image_contrast < mean_contrast - std_contrast:
        target_std = mean_contrast - std_contrast
        image = adjust_contrast(image, target_std)
    elif image_contrast > mean_contrast + std_contrast:
        target_std = mean_contrast + std_contrast
        image = adjust_contrast(image, target_std)

    normalized_images.append(image)

for image_name, image in zip(file_names, normalized_images):
    cv2.imwrite(image_name, image)

