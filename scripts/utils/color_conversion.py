import pandas as pd
import matplotlib

# Read the CSV into a DataFrame
colors_df = pd.read_csv("scripts/utils/colors.csv")

# Create dictionaries for fast look-up
color_to_rgb = dict(zip(colors_df["ColorName"].str.lower(), colors_df["RGB"]))
color_to_hex = dict(zip(colors_df["ColorName"].str.lower(), colors_df["HexColor"]))


def get_rgb(color_name):
    color_name = color_name.lower()
    if color_name in color_to_rgb:
        return color_to_rgb[color_name]
    try:
        # Convert using matplotlib and then return the RGB values as a string
        rgba = matplotlib.colors.to_rgba(color_name)
        return ",".join(map(str, rgba[:-1]))
    except ValueError:
        return None


def get_hex(color_name):
    color_name = color_name.lower()
    if color_name in color_to_hex:
        return color_to_hex[color_name]
    try:
        # Convert using matplotlib
        return matplotlib.colors.to_hex(color_name)
    except ValueError:
        return None


# # Examples
# print(get_rgb("darkgrey"))  # Should give RGB values
# print(get_hex("darkgrey"))  # Should give HEX values
# print(get_rgb("invalidcolor"))  # Should return None
# print(get_hex("invalidcolor"))  # Should return None
