import matplotlib.pyplot as plt
from PIL import Image

img1 = Image.open('TaC_SNAP/GsfCurves/GSF_both.png')
img2 = Image.open('TaC_ACE/GsfCurves/GSF_both.png')
img3 = Image.open('TaC_MEAM/GsfCurves/GSF_both.png')

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

axes[0].imshow(img1)
axes[0].axis('off')
axes[0].set_title('SNAP')

axes[1].imshow(img2)
axes[1].axis('off')
axes[1].set_title('ACE')

axes[2].imshow(img3)
axes[2].axis('off')
axes[2].set_title('MEAM')

plt.tight_layout()
plt.savefig('GSF_comparison.png', dpi=300)
plt.show()
