import os
import cv2
import numpy as np
import pandas as pd
from skimage import color
import pywt
from scipy import stats


def gaussian_pyramid(image, levels=3):
    """Create a Gaussian pyramid of the image."""
    pyramid = [image]
    for _ in range(1, levels):
        image = cv2.pyrDown(image)
        pyramid.append(image)
    return pyramid


def compute_contrast_energy(luminance, sigma1=1.0, sigma2=2.0):
    """Compute contrast energy using a difference of Gaussians."""
    gaussian1 = cv2.GaussianBlur(luminance, (0, 0), sigmaX=sigma1)
    gaussian2 = cv2.GaussianBlur(luminance, (0, 0), sigmaX=sigma2)
    contrast_energy = (gaussian1 - gaussian2) ** 2
    return contrast_energy


def compute_orientation_energy(luminance, ksize=3):
    """Compute orientation energy using oriented opponent energy."""
    gradients = []
    for theta in np.linspace(0, np.pi, 4, endpoint=False):
        kernel = cv2.getGaborKernel((ksize, ksize), sigma=1.0, theta=theta, lambd=10.0, gamma=0.5)
        filtered = cv2.filter2D(luminance, cv2.CV_64F, kernel)
        gradients.append(filtered)
    gradients = np.array(gradients)
    orientation_energy = np.sum(gradients ** 2, axis=0)
    return orientation_energy


def compute_local_covariance(feature_maps):
    """Compute local covariance of the feature maps."""
    covariance_maps = []
    for feature_map in feature_maps:
        mean_map = cv2.GaussianBlur(feature_map, (0, 0), sigmaX=1.0)
        variance_map = cv2.GaussianBlur((feature_map - mean_map) ** 2, (0, 0), sigmaX=1.0)
        covariance_maps.append(np.sqrt(variance_map))
    return covariance_maps


def resize_to_max_shape(feature_maps, target_shape):
    """Resize feature maps to the target shape."""
    resized_maps = [cv2.resize(fm, (target_shape[1], target_shape[0])) for fm in feature_maps]
    return resized_maps


def combine_scales(feature_pyramids, target_shape):
    """Combine feature maps across scales by taking the maximum at each pixel."""
    combined_maps = []
    for feature_maps in zip(*feature_pyramids):
        resized_maps = resize_to_max_shape(feature_maps, target_shape)
        combined_map = np.max(np.stack(resized_maps, axis=-1), axis=-1)
        combined_maps.append(combined_map)
    return combined_maps


def normalize_feature_maps(feature_maps, epsilon=1e-10):
    """Normalize feature maps using the standard deviation across a range of images."""
    normalized_maps = []
    for feature_map in feature_maps:
        std_dev = np.std(feature_map)
        normalized_map = feature_map / (std_dev + epsilon)
        normalized_maps.append(normalized_map)
    return normalized_maps


def compute_clutter_score(color_clutter, contrast_clutter, orientation_clutter):
    """Combine clutter maps into a single clutter score."""
    combined_clutter = np.cbrt(color_clutter) + np.sqrt(contrast_clutter) + np.sqrt(orientation_clutter)
    return np.mean(combined_clutter)


def feature_congestion_clutter(image):
    """Compute the Feature Congestion measure of visual clutter."""

    # Convert the image to the CIELab color space
    lab_image = color.rgb2lab(image)

    # Extract the L (luminance) channel for processing
    luminance = lab_image[:, :, 0]

    # Create a Gaussian pyramid
    pyramid = gaussian_pyramid(lab_image)

    # Extract features at each scale
    color_pyramids = [p[:, :, 1:] for p in pyramid]  # Extract a*b* channels
    contrast_pyramids = [compute_contrast_energy(p[:, :, 0]) for p in pyramid]
    orientation_pyramids = [compute_orientation_energy(p[:, :, 0]) for p in pyramid]

    # Determine the target shape (original image shape)
    target_shape = pyramid[0].shape[:2]

    # Compute local covariance (clutter) for each feature
    color_clutter = compute_local_covariance(color_pyramids)
    contrast_clutter = compute_local_covariance(contrast_pyramids)
    orientation_clutter = compute_local_covariance(orientation_pyramids)

    # Combine across scales
    combined_color = combine_scales(color_clutter, target_shape)
    combined_contrast = combine_scales(contrast_clutter, target_shape)
    combined_orientation = combine_scales(orientation_clutter, target_shape)

    # Normalize the feature maps
    normalized_color = normalize_feature_maps(combined_color)
    normalized_contrast = normalize_feature_maps(combined_contrast)
    normalized_orientation = normalize_feature_maps(combined_orientation)

    # Compute the final clutter score
    clutter_score = compute_clutter_score(normalized_color, normalized_contrast, normalized_orientation)

    return clutter_score


def compute_shannon_entropy(coefficients, num_bins=256):
    """Compute the Shannon entropy of wavelet coefficients within a subband."""
    histogram, _ = np.histogram(coefficients, bins=num_bins, density=True)
    histogram = histogram[histogram > 0]  # Exclude zero probabilities
    entropy = -np.sum(histogram * np.log2(histogram))
    return entropy


def wavelet_decomposition(channel, wavelet='db1', level=5):
    """Decompose a channel into wavelet subbands using PyWavelets."""
    coeffs = pywt.wavedec2(channel, wavelet=wavelet, level=level)
    return coeffs


def compute_subband_entropy(image):
    """Compute the Subband Entropy measure of visual clutter."""

    # Convert the image to the CIELab color space
    lab_image = color.rgb2lab(image)

    # Extract the L (luminance) and a, b (chrominance) channels
    L_channel = lab_image[:, :, 0]
    a_channel = lab_image[:, :, 1]
    b_channel = lab_image[:, :, 2]

    # Decompose each channel into wavelet subbands
    L_coeffs = wavelet_decomposition(L_channel)
    a_coeffs = wavelet_decomposition(a_channel)
    b_coeffs = wavelet_decomposition(b_channel)

    # Compute the entropy for each subband (excluding the approximation coefficients)
    L_entropy = sum(compute_shannon_entropy(coeff) for coeff in L_coeffs[1:])
    a_entropy = sum(compute_shannon_entropy(coeff) for coeff in a_coeffs[1:])
    b_entropy = sum(compute_shannon_entropy(coeff) for coeff in b_coeffs[1:])

    # Compute the weighted sum of luminance and chrominance entropies
    subband_entropy = 0.84 * L_entropy + 0.08 * a_entropy + 0.08 * b_entropy

    return subband_entropy


def compute_luminance(image):
    """
    Compute the overall luminance of an image using the CIELab color space.

    Parameters:
        image (numpy array): Input image in RGB format.

    Returns:
        float: Overall luminance of the image.
    """
    # Convert the image to CIELab color space
    lab_image = color.rgb2lab(image)

    # Extract the L channel (luminance)
    luminance = lab_image[:, :, 0]

    # Calculate the mean luminance
    mean_luminance = np.mean(luminance)

    return mean_luminance


def compute_contrast(image):
    """
    Compute the contrast of an image using the standard deviation of luminance in CIELab color space (similar to GIMP).

    Parameters:
        image (numpy array): Input image in RGB format.

    Returns:
        float: Contrast of the image.
    """
    # Convert the image to CIELab color space
    lab_image = color.rgb2lab(image)

    # Extract the L channel (luminance)
    luminance = lab_image[:, :, 0]

    # Calculate the contrast (standard deviation of luminance)
    contrast = np.std(luminance)

    return contrast


def add_gray_background(image, background_color=(128, 128, 128)):
    """
    Add a gray background to an image with transparency.

    Parameters:
        image (numpy array): Input image with RGBA channels.
        background_color (tuple): Background color in RGB format (default is gray #808080).

    Returns:
        numpy array: Image with the gray background added.
    """
    # Check if the image has an alpha channel
    if image.shape[2] == 4:
        # Separate the color and alpha channels
        bgr = image[:, :, :3]
        alpha = image[:, :, 3] / 255.0

        # Create a background image with the specified background color
        background = np.full_like(bgr, background_color, dtype=np.uint8)

        # Blend the image with the background
        image_with_background = (alpha[..., None] * bgr + (1 - alpha[..., None]) * background).astype(np.uint8)

        return image_with_background
    else:
        # If no alpha channel, return the original image
        return image


# Compute complexity scores for stimuli of experiment 1:
print("Stimuli of Experiment 1")
scores_exp1 = pd.DataFrame()
path = os.path.join('Study 1', 'Experiment', 'gaze_avoidance_task', 'stimuli')
files = [file for file in os.listdir(path) if file.endswith(".png") and not "Kopie" in file]
for image_file in files:
    # image_file = files[0]
    image = cv2.imread(os.path.join(path, image_file), cv2.IMREAD_UNCHANGED)

    # Add gray background
    image = add_gray_background(image)

    scores_exp1 = pd.concat([scores_exp1, pd.DataFrame({"file": [image_file],
                                              "clutter_score": [feature_congestion_clutter(image)],
                                              "entropy_score": [compute_subband_entropy(image)],
                                              "luminance": [compute_luminance(image)],
                                              "contrast": [compute_contrast(image)]})])

t_clutter, p_clutter = stats.ttest_ind(scores_exp1.loc[scores_exp1["file"].str.contains("Nonsocial"), "clutter_score"], scores_exp1.loc[~scores_exp1["file"].str.contains("Nonsocial"), "clutter_score"])
print(f"Clutter: t = {t_clutter}, p = {p_clutter}")
t_entropy, p_entropy = stats.ttest_ind(scores_exp1.loc[scores_exp1["file"].str.contains("Nonsocial"), "entropy_score"], scores_exp1.loc[~scores_exp1["file"].str.contains("Nonsocial"), "entropy_score"])
print(f"Entropy: t = {t_entropy}, p = {p_entropy}")
t_luminance, p_luminance = stats.ttest_ind(scores_exp1.loc[scores_exp1["file"].str.contains("Nonsocial"), "luminance"], scores_exp1.loc[~scores_exp1["file"].str.contains("Nonsocial"), "luminance"])
print(f"Luminance: t = {t_luminance}, p = {p_luminance}")
t_contrast, p_contrast = stats.ttest_ind(scores_exp1.loc[scores_exp1["file"].str.contains("Nonsocial"), "contrast"], scores_exp1.loc[~scores_exp1["file"].str.contains("Nonsocial"), "contrast"])
print(f"Contrast: t = {t_contrast}, p = {p_contrast}")

# Compute complexity scores for stimuli of experiment 2:
print("Stimuli of Experiment 2")
scores_exp2 = pd.DataFrame()
path = os.path.join('Study 2', 'Experiment', 'attentional_competition_task', 'stimuli')
files = [file for file in os.listdir(path) if file.endswith(".png") and not "Kopie" in file]
for image_file in files:
    # image_file = files[0]
    image = cv2.imread(os.path.join(path, image_file), cv2.IMREAD_UNCHANGED)

    # Add gray background
    image = add_gray_background(image)

    scores_exp2 = pd.concat([scores_exp2, pd.DataFrame({"file": [image_file],
                                              "clutter_score": [feature_congestion_clutter(image)],
                                              "entropy_score": [compute_subband_entropy(image)],
                                              "luminance": [compute_luminance(image)],
                                              "contrast": [compute_contrast(image)]})])

t_clutter, p_clutter = stats.ttest_ind(scores_exp2.loc[scores_exp2["file"].str.contains("DALLE"), "clutter_score"], scores_exp2.loc[~scores_exp2["file"].str.contains("DALLE"), "clutter_score"])
print(f"Clutter: t = {t_clutter}, p = {p_clutter}")
t_entropy, p_entropy = stats.ttest_ind(scores_exp2.loc[scores_exp2["file"].str.contains("DALLE"), "entropy_score"], scores_exp2.loc[~scores_exp2["file"].str.contains("DALLE"), "entropy_score"])
print(f"Entropy: t = {t_entropy}, p = {p_entropy}")
t_luminance, p_luminance = stats.ttest_ind(scores_exp2.loc[scores_exp2["file"].str.contains("DALLE"), "luminance"], scores_exp2.loc[~scores_exp2["file"].str.contains("DALLE"), "luminance"])
print(f"Luminance: t = {t_luminance}, p = {p_luminance}")
t_contrast, p_contrast = stats.ttest_ind(scores_exp2.loc[scores_exp2["file"].str.contains("DALLE"), "contrast"], scores_exp2.loc[~scores_exp2["file"].str.contains("DALLE"), "contrast"])
print(f"Contrast: t = {t_contrast}, p = {p_contrast}")
