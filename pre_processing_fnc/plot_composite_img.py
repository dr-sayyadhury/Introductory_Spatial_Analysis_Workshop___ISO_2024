import numpy as np
import matplotlib.pyplot as plt


def plot_composite_image(image, colors=None):
    """
    Plot a composite image with the given colors for each channel.
    :param image: The composite image to plot.
    :param colors: A list of colors for each channel.
    """
 
    # Check the data's shape and maximum value
    print("Shape of iF_crop:", image.shape)
    for i in range(4):
        print(f"Channel {i} max: {np.max(image[i])}, min: {np.min(image[i])}")

    if colors is None:
        # Define colors, with DAPI visualized in blue
        colors = [
            np.array([0, 0, 1]),  # Blue for DAPI
            np.array([1, 0, 0]),  # Red
            np.array([0, 1, 0]),  # Green
            np.array([1, 1, 0])   # Yellow for the fourth channel
        ]
    else:
        # Ensure the colors are numpy arrays
        colors = [np.array(color) for color in colors]

    # Create a composite image with 3 channels (RGB)
    composite_image = np.zeros((image.shape[1], image.shape[2], 3), dtype=np.float32)

    # Assign each channel to a color in RGB, using a more conservative normalization approach
    for i in range(4):
        channel_data = image[i, :, :]
        if np.max(channel_data) > 0:  # Avoid division by zero and unnecessary normalization
            normalized_data = channel_data / np.max(channel_data)
            composite_image += normalized_data[:, :, np.newaxis] * colors[i]

    return composite_image
#composite_image += channel_data[:, :, np.newaxis] * colors[i]


# Ensure the image is normalized to the maximum of the composite_image
#if np.max(composite_image) > 0:
#    composite_image /= np.max(composite_image)
