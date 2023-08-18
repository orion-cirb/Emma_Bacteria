# Emma_Bacteria

* **Developed for:** Emma
* **Team:** Espeli
* **Date:** July 2023
* **Software:** Fiji


### Images description

3D images taken with a x60 objective

2 channels:
  1. *488:* Gene of interest
  2. *561:* mCherry-marked bacteria

### Plugin description

* Perform average intensity Z-projection of both channels
* Detect bacteria on channel 2 with Omnipose
* Measure integrated intensity of each bacterium in channels 2 and 3, and give their ratio

### Dependencies

* **3DImageSuite** Fiji plugin
* **Omnipose** conda environment + *bact_fluor_omnitorch_0* model

### Version history

Version 1 released on July 28, 2023.

