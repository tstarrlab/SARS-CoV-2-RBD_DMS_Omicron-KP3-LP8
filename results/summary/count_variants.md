# Count variants in each sample
This Python Jupyter notebook counts occurrences of each barcode in each sample from Illumina barcode sequencing, and adds these counts to the codon variant table.

## Set up analysis
### Import Python modules.
Use [plotnine](https://plotnine.readthedocs.io/en/stable/) for ggplot2-like plotting.

The analysis relies heavily on the Bloom lab's [dms_variants](https://jbloomlab.github.io/dms_variants) package:


```python
import itertools
import multiprocessing
import multiprocessing.pool
import os
import warnings

import alignparse
import alignparse.targets

import dms_variants.codonvarianttable
from dms_variants.constants import CBPALETTE
import dms_variants.illuminabarcodeparser
import dms_variants.utils
import dms_variants.plotnine_themes

from IPython.display import display, HTML

import pandas as pd

from plotnine import *

import yaml
```

Set [plotnine](https://plotnine.readthedocs.io/en/stable/) theme to the gray-grid one defined in `dms_variants`:


```python
theme_set(dms_variants.plotnine_themes.theme_graygrid())
```

Versions of key software:


```python
print(f"Using alignparse version {alignparse.__version__}")
print(f"Using dms_variants version {dms_variants.__version__}")
```

    Using alignparse version 0.6.0
    Using dms_variants version 1.4.3


Ignore warnings that clutter output:


```python
warnings.simplefilter('ignore')
```

### Parameters for notebook
Read the configuration file:


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)
```

Make output directory if needed:


```python
os.makedirs(config['counts_dir'], exist_ok=True)
os.makedirs(config['figs_dir'], exist_ok=True)
```

## Input variant tables
Initialize the table of barcode-variant pairs from the respective `process_ccs` notebooks for each background.


```python
variants = pd.read_csv(config['codon_variant_table_file_KP3'], na_filter=None)
variants = variants.append(pd.read_csv(config['codon_variant_table_file_LP8'], na_filter=None))

variants = variants.reset_index(drop=True)

display(HTML(variants.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>target</th>
      <th>library</th>
      <th>barcode</th>
      <th>variant_call_support</th>
      <th>codon_substitutions</th>
      <th>aa_substitutions</th>
      <th>n_codon_substitutions</th>
      <th>n_aa_substitutions</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>KP3</td>
      <td>pool1</td>
      <td>AAAAAAAAAAAGAAAA</td>
      <td>1</td>
      <td>AAG148ATT</td>
      <td>K148I</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>KP3</td>
      <td>pool1</td>
      <td>AAAAAAAAAAGGAAAC</td>
      <td>1</td>
      <td>AAG110ACT</td>
      <td>K110T</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>KP3</td>
      <td>pool1</td>
      <td>AAAAAAAAAGAGGGAT</td>
      <td>21</td>
      <td>AAC107CCA</td>
      <td>N107P</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>KP3</td>
      <td>pool1</td>
      <td>AAAAAAAAAGAGTACA</td>
      <td>6</td>
      <td>ATT138ATG</td>
      <td>I138M</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>KP3</td>
      <td>pool1</td>
      <td>AAAAAAAAATACCAGT</td>
      <td>4</td>
      <td>ACC3GGT</td>
      <td>T3G</td>
      <td>1</td>
      <td>1</td>
    </tr>
  </tbody>
</table>


Are there any barcodes in the same library that are shared across targets?
If so, we need to get rid of those as they will be confounded in barcode parsing:


```python
dup_barcodes = (
    variants
    .groupby(['library', 'barcode'])
    .size()
    .rename('duplicate_count')
    .reset_index()
    .query('duplicate_count > 1')
    )

print('Here are duplicated barcodes:')
display(HTML(dup_barcodes.head().to_html(index=False)))

print(f"\nRemoving the {len(dup_barcodes)} duplicated barcodes."
      f"Started with {len(variants)} barcodes:")
variants = (
    variants
    .merge(dup_barcodes, on=['library', 'barcode'], how='outer')
    .query('duplicate_count.isnull()', engine='python')
    )
print(f"After removing duplicates, there are {len(variants)} barcodes.")
```

    Here are duplicated barcodes:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>barcode</th>
      <th>duplicate_count</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>pool1</td>
      <td>AAACCCAAAAATGATA</td>
      <td>2</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>AAGCAATATAGACATT</td>
      <td>2</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>ACCAAATTTGGCAATA</td>
      <td>2</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>ACTGAAAAATAGAGAA</td>
      <td>2</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>ATCATAAATAACCAAT</td>
      <td>2</td>
    </tr>
  </tbody>
</table>


    
    Removing the 24 duplicated barcodes.Started with 408175 barcodes:
    After removing duplicates, there are 408127 barcodes.


Pull out a target sequence for matching to the barcode and flanking sequence regions. Note, in this pipeline this is ok because our different backgrounds don't have differing flanks or other features within the actual N16 region covered in Illumina sequencing. If ever placing in-line barcodes here in the future, we would need to modify this.


```python
# get wildtype gene sequence for primary target
targets = alignparse.targets.Targets(seqsfile=config['amplicons_KP3'],
                                     feature_parse_specs=config['feature_parse_specs_KP3'])
```

## Setup to parse barcodes
Read data frame with list of all barcode runs.


```python
# barcode runs with R1 files by semicolon string split
barcode_runs = (pd.read_csv(config['barcode_runs'])
                .assign(R1=lambda x: x['R1'].str.split('; '))
                )

display(HTML(barcode_runs.to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>sample</th>
      <th>sample_type</th>
      <th>sort_bin</th>
      <th>concentration</th>
      <th>date</th>
      <th>number_cells</th>
      <th>R1</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>1.0</td>
      <td>250626</td>
      <td>2069353.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250626-exp1-s1-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_01_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>1.0</td>
      <td>250626</td>
      <td>1082207.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250626-exp1-s1-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_01_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>1.0</td>
      <td>250626</td>
      <td>2140833.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250626-exp1-s1-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_01_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>1.0</td>
      <td>250626</td>
      <td>2318187.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250626-exp1-s1-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_02_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>2.0</td>
      <td>250626</td>
      <td>1970977.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250626-exp1-s2-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_02_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>2.0</td>
      <td>250626</td>
      <td>1153905.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250626-exp1-s2-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_02_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>2.0</td>
      <td>250626</td>
      <td>2231043.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250626-exp1-s2-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_02_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>2.0</td>
      <td>250626</td>
      <td>1361470.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250626-exp1-s2-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_03_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>3.0</td>
      <td>250626</td>
      <td>2411129.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250626-exp1-s3-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_03_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>3.0</td>
      <td>250626</td>
      <td>1365576.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250626-exp1-s3-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_03_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>3.0</td>
      <td>250626</td>
      <td>2207136.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250626-exp1-s3-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_03_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>3.0</td>
      <td>250626</td>
      <td>581304.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250626-exp1-s3-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_04_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>4.0</td>
      <td>250626</td>
      <td>3982304.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250626-exp1-s4-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_04_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>4.0</td>
      <td>250626</td>
      <td>2046138.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250626-exp1-s4-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_04_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>4.0</td>
      <td>250626</td>
      <td>588604.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250626-exp1-s4-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_04_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>4.0</td>
      <td>250626</td>
      <td>10089.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250626-exp1-s4-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_05_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>5.0</td>
      <td>250626</td>
      <td>5737668.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250626-exp1-s5-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_05_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>5.0</td>
      <td>250626</td>
      <td>1333844.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250626-exp1-s5-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_05_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>5.0</td>
      <td>250626</td>
      <td>64235.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250626-exp1-s5-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_05_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>5.0</td>
      <td>250626</td>
      <td>1612.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250626-exp1-s5-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_06_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>6.0</td>
      <td>250626</td>
      <td>6615270.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250626-exp1-s6-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_06_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>6.0</td>
      <td>250626</td>
      <td>216074.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250626-exp1-s6-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_06_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>6.0</td>
      <td>250626</td>
      <td>5535.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250626-exp1-s6-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_06_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>6.0</td>
      <td>250626</td>
      <td>1871.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250626-exp1-s6-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_07_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>7.0</td>
      <td>250626</td>
      <td>6813335.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250626-exp1-s7-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_07_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>7.0</td>
      <td>250626</td>
      <td>198931.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250626-exp1-s7-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_07_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>7.0</td>
      <td>250626</td>
      <td>1355.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250626-exp1-s7-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_07_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>7.0</td>
      <td>250626</td>
      <td>1159.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250626-exp1-s7-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_08_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>8.0</td>
      <td>250626</td>
      <td>6849588.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250626-exp1-s8-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_08_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>8.0</td>
      <td>250626</td>
      <td>197657.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250626-exp1-s8-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_08_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>8.0</td>
      <td>250626</td>
      <td>1727.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250626-exp1-s8-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_08_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>8.0</td>
      <td>250626</td>
      <td>1247.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250626-exp1-s8-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_09_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>9.0</td>
      <td>250626</td>
      <td>6510069.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250626-exp1-s9-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_09_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>9.0</td>
      <td>250626</td>
      <td>181654.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250626-exp1-s9-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_09_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>9.0</td>
      <td>250626</td>
      <td>2174.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250626-exp1-s9-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>TiteSeq_09_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>9.0</td>
      <td>250626</td>
      <td>1590.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250626-exp1-s9-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_01_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>1.0</td>
      <td>250703</td>
      <td>1386204.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250703-exp2-s1-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_01_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>1.0</td>
      <td>250703</td>
      <td>785808.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250703-exp2-s1-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_01_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>1.0</td>
      <td>250703</td>
      <td>1962408.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250703-exp2-s1-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_01_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>1.0</td>
      <td>250703</td>
      <td>2283804.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250703-exp2-s1-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_02_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>2.0</td>
      <td>250703</td>
      <td>1555783.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250703-exp2-s2-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_02_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>2.0</td>
      <td>250703</td>
      <td>1121256.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250703-exp2-s2-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_02_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>2.0</td>
      <td>250703</td>
      <td>3073945.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250703-exp2-s2-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_02_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>2.0</td>
      <td>250703</td>
      <td>937363.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250703-exp2-s2-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_03_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>3.0</td>
      <td>250703</td>
      <td>1919464.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250703-exp2-s3-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_03_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>3.0</td>
      <td>250703</td>
      <td>1165430.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250703-exp2-s3-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_03_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>3.0</td>
      <td>250703</td>
      <td>2618732.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250703-exp2-s3-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_03_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>3.0</td>
      <td>250703</td>
      <td>966612.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250703-exp2-s3-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_04_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>4.0</td>
      <td>250703</td>
      <td>2473369.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250703-exp2-s4-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_04_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>4.0</td>
      <td>250703</td>
      <td>1898680.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250703-exp2-s4-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_04_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>4.0</td>
      <td>250703</td>
      <td>2026836.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250703-exp2-s4-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_04_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>4.0</td>
      <td>250703</td>
      <td>169464.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250703-exp2-s4-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_05_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>5.0</td>
      <td>250703</td>
      <td>3955061.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250703-exp2-s5-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_05_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>5.0</td>
      <td>250703</td>
      <td>1920276.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250703-exp2-s5-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_05_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>5.0</td>
      <td>250703</td>
      <td>383892.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250703-exp2-s5-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_05_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>5.0</td>
      <td>250703</td>
      <td>31493.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250703-exp2-s5-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_06_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>6.0</td>
      <td>250703</td>
      <td>4663633.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250703-exp2-s6-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_06_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>6.0</td>
      <td>250703</td>
      <td>288008.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250703-exp2-s6-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_06_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>6.0</td>
      <td>250703</td>
      <td>8190.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250703-exp2-s6-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_06_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>6.0</td>
      <td>250703</td>
      <td>10063.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250703-exp2-s6-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_07_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>7.0</td>
      <td>250703</td>
      <td>4828432.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250703-exp2-s7-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_07_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>7.0</td>
      <td>250703</td>
      <td>379187.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250703-exp2-s7-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_07_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>7.0</td>
      <td>250703</td>
      <td>14327.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250703-exp2-s7-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_07_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>7.0</td>
      <td>250703</td>
      <td>14691.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250703-exp2-s7-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_08_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>8.0</td>
      <td>250703</td>
      <td>6067314.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250703-exp2-s8-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_08_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>8.0</td>
      <td>250703</td>
      <td>402007.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250703-exp2-s8-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_08_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>8.0</td>
      <td>250703</td>
      <td>17314.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250703-exp2-s8-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_08_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>8.0</td>
      <td>250703</td>
      <td>17209.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250703-exp2-s8-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_09_bin1</td>
      <td>TiteSeq</td>
      <td>1</td>
      <td>9.0</td>
      <td>250703</td>
      <td>6534658.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250703-exp2-s9-b1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_09_bin2</td>
      <td>TiteSeq</td>
      <td>2</td>
      <td>9.0</td>
      <td>250703</td>
      <td>188870.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250703-exp2-s9-b2_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_09_bin3</td>
      <td>TiteSeq</td>
      <td>3</td>
      <td>9.0</td>
      <td>250703</td>
      <td>914.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250703-exp2-s9-b3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>TiteSeq_09_bin4</td>
      <td>TiteSeq</td>
      <td>4</td>
      <td>9.0</td>
      <td>250703</td>
      <td>907.000</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250703-exp2-s9-b4_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>SortSeq_bin1</td>
      <td>SortSeq</td>
      <td>1</td>
      <td>NaN</td>
      <td>250731</td>
      <td>3293761.680</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250731-exp3rep1-RBD-b1-1_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250731-exp3rep1-RBD-b1-2_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250731-exp3rep1-RBD-b1-3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>SortSeq_bin2</td>
      <td>SortSeq</td>
      <td>2</td>
      <td>NaN</td>
      <td>250731</td>
      <td>2622277.310</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250731-exp3rep1-RBD-b2-1_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250731-exp3rep1-RBD-b2-2_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250731-exp3rep1-RBD-b2-3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>SortSeq_bin3</td>
      <td>SortSeq</td>
      <td>3</td>
      <td>NaN</td>
      <td>250731</td>
      <td>2228095.320</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250731-exp3rep1-RBD-b3-1_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250731-exp3rep1-RBD-b3-2_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250731-exp3rep1-RBD-b3-3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>SortSeq_bin4</td>
      <td>SortSeq</td>
      <td>4</td>
      <td>NaN</td>
      <td>250731</td>
      <td>2553606.720</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250731-exp3rep1-RBD-b4-1_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250731-exp3rep1-RBD-b4-2_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250731-exp3rep1-RBD-b4-3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>SortSeq_bin1</td>
      <td>SortSeq</td>
      <td>1</td>
      <td>NaN</td>
      <td>250731</td>
      <td>3197897.083</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250731-exp3rep2-RBD-b1-1_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250731-exp3rep2-RBD-b1-2_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250731-exp3rep2-RBD-b1-3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>SortSeq_bin2</td>
      <td>SortSeq</td>
      <td>2</td>
      <td>NaN</td>
      <td>250731</td>
      <td>3421500.160</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250731-exp3rep2-RBD-b2-1_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250731-exp3rep2-RBD-b2-2_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250731-exp3rep2-RBD-b2-3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>SortSeq_bin3</td>
      <td>SortSeq</td>
      <td>3</td>
      <td>NaN</td>
      <td>250731</td>
      <td>4067189.555</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250731-exp3rep2-RBD-b3-1_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250731-exp3rep2-RBD-b3-2_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250731-exp3rep2-RBD-b3-3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>pool2</td>
      <td>SortSeq_bin4</td>
      <td>SortSeq</td>
      <td>4</td>
      <td>NaN</td>
      <td>250731</td>
      <td>4103152.695</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250731-exp3rep2-RBD-b4-1_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250731-exp3rep2-RBD-b4-2_R1_001.fastq.gz, /uufs/chpc.utah.edu/common/home/starr-group1/sequencing/ALT/2025/250821_illumina_SARS2-kp3-lp8/250731-exp3rep2-RBD-b4-3_R1_001.fastq.gz]</td>
    </tr>
  </tbody>
</table>


Make sure library / sample combinations are unique:


```python
assert len(barcode_runs) == len(barcode_runs.groupby(['library', 'sample']))
```

Make sure the the libraries for which we have barcode runs are all in our variant table:


```python
unknown_libs = set(barcode_runs['library']) - set(variants['library'])
if unknown_libs:
    raise ValueError(f"Libraries with barcode runs not in variant table: {unknown_libs}")
```

Now we initialize an [IlluminaBarcodeParser](https://jbloomlab.github.io/dms_variants/dms_variants.illuminabarcodeparser.html#dms_variants.illuminabarcodeparser.IlluminaBarcodeParser) for each library.

First, get the length of the barcode from the alignment target after making sure the same length for all targets:


```python
bclen = len(targets.targets[0].get_feature('barcode').seq)

assert (bclen == len(target.get_feature('barcode').seq) for target in targets.targets)

print(f"Barcodes of length {bclen}")
```

    Barcodes of length 16


The other barcode parsing params come from the config file:


```python
parser_params = config['illumina_barcode_parser_params']

display(HTML(
    pd.Series(parser_params, name='value')
    .rename_axis(index='parameter')
    .reset_index()
    .to_html(index=False)
    ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>parameter</th>
      <th>value</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>upstream</td>
      <td>GGCCGC</td>
    </tr>
    <tr>
      <td>downstream</td>
      <td></td>
    </tr>
    <tr>
      <td>minq</td>
      <td>20</td>
    </tr>
    <tr>
      <td>upstream_mismatch</td>
      <td>1</td>
    </tr>
    <tr>
      <td>downstream_mismatch</td>
      <td>0</td>
    </tr>
  </tbody>
</table>


The parser needs to know the set of valid barcodes, which are stored in the variant table and are different for each library.
So we create a different parser for each library using these valid barcode sets:


```python
# create dict keyed by library, value is parser for library
parsers = {lib: dms_variants.illuminabarcodeparser.IlluminaBarcodeParser(
                    bclen=bclen,
                    valid_barcodes=variants.loc[variants['library']==lib]['barcode'],
                    **parser_params)
           for lib in set(variants['library'])}

print('Number of valid barcodes searched for by each parser:')
display(HTML(
    pd.DataFrame([(lib, len(p.valid_barcodes)) for lib, p in parsers.items()],
                 columns=['library', 'number of valid barcodes'])
    .to_html(index=False)
    ))
```

    Number of valid barcodes searched for by each parser:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>number of valid barcodes</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>pool2</td>
      <td>203512</td>
    </tr>
    <tr>
      <td>pool1</td>
      <td>204615</td>
    </tr>
  </tbody>
</table>


## Parse barcodes
We now parse the barcodes.
Since this will take a while, we utilize multiple CPUs via the Python [multiprocessing](https://docs.python.org/3.6/library/multiprocessing.html) module.
First, determine how many CPUs to use.
We use the minimum of the user-specified number hardcoded below and the number actually available.
(If you are running *interactively* on the Hutch cluster, you may need to reduce the number below in order to avoid an error as there is an enforced CPU limit on the home `rhino` nodes):


```python
ncpus = min(config['max_cpus'], multiprocessing.cpu_count())
print(f"Using {ncpus} CPUs")
```

    Using 8 CPUs


Parse the barcodes in parallel via a [multiprocessing.Pool](https://docs.python.org/3.6/library/multiprocessing.html#multiprocessing.pool.Pool) using all the available CPUs to get a list of the data frames with barcode counts / fates for each sample:


```python
def process_func(parser, r1files, library, sample):
    """Convenience function to be starmapped to multiprocessing pool."""
    return parser.parse(r1files, add_cols={'library': library, 'sample': sample})

# parallel computation of list of data frames
with multiprocessing.pool.Pool(processes=ncpus) as pool:
    bclist = pool.starmap(
                process_func,
                [(parsers[run.library], run.R1, run.library, run.sample)
                  for run in barcode_runs.itertuples()],
                )
```

Now concatenate the list into data frames of barcode counts and barcode fates:


```python
counts = pd.concat([samplecounts for samplecounts, _ in bclist],
                   sort=False,
                   ignore_index=True)

print('First few lines of counts data frame:')
display(HTML(counts.head().to_html(index=False)))

fates = pd.concat([samplefates for _, samplefates in bclist],
                  sort=False,
                  ignore_index=True)

print('First few lines of fates data frame:')
display(HTML(fates.head().to_html(index=False)))
```

    First few lines of counts data frame:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>barcode</th>
      <th>count</th>
      <th>library</th>
      <th>sample</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>TACCGAACCGGTGCAA</td>
      <td>2861</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>ACCCACACAATGCAGG</td>
      <td>2659</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>ATTACGTTTTTACCTC</td>
      <td>2571</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>TAACATATGGAACACG</td>
      <td>2320</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>CTTCGGAACGCTAAAC</td>
      <td>2319</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
  </tbody>
</table>


    First few lines of fates data frame:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>fate</th>
      <th>count</th>
      <th>library</th>
      <th>sample</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>valid barcode</td>
      <td>9171735</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>invalid barcode</td>
      <td>1553051</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>low quality barcode</td>
      <td>801500</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>unparseable barcode</td>
      <td>197382</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
    <tr>
      <td>failed chastity filter</td>
      <td>0</td>
      <td>pool1</td>
      <td>TiteSeq_01_bin1</td>
    </tr>
  </tbody>
</table>


## Examine fates of parsed barcodes
First, we'll analyze the "fates" of the parsed barcodes.
These fates represent what happened to each Illumina read we parsed:
 - Did the barcode read fail the Illumina chastity filter?
 - Was the barcode *unparseable* (i.e., the read didn't appear to be a valid barcode based on flanking regions)?
 - Was the barcode sequence too *low quality* based on the Illumina quality scores?
 - Was the barcode parseable but *invalid* (i.e., not in our list of variant-associated barcodes in the codon variant table)?
 - Was the barcode *valid*, and so will be added to variant counts.
 
First, we just write a CSV file with all the barcode fates:


```python
fatesfile = os.path.join(config['counts_dir'], 'barcode_fates.csv')
print(f"Writing barcode fates to {fatesfile}")
fates.to_csv(fatesfile, index=False)
```

    Writing barcode fates to results/counts/barcode_fates.csv


Next, we tabulate the barcode fates in wide format:


```python
display(HTML(fates
             .pivot_table(columns='fate',
                          values='count',
                          index=['library', 'sample'])
             .to_html()
             ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>fate</th>
      <th>failed chastity filter</th>
      <th>invalid barcode</th>
      <th>low quality barcode</th>
      <th>unparseable barcode</th>
      <th>valid barcode</th>
    </tr>
    <tr>
      <th>library</th>
      <th>sample</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="40" valign="top">pool1</th>
      <th>SortSeq_bin1</th>
      <td>0</td>
      <td>2640678</td>
      <td>1286253</td>
      <td>354848</td>
      <td>13484014</td>
    </tr>
    <tr>
      <th>SortSeq_bin2</th>
      <td>0</td>
      <td>2118462</td>
      <td>1201927</td>
      <td>267259</td>
      <td>12324497</td>
    </tr>
    <tr>
      <th>SortSeq_bin3</th>
      <td>0</td>
      <td>2094050</td>
      <td>1267751</td>
      <td>267374</td>
      <td>12409077</td>
    </tr>
    <tr>
      <th>SortSeq_bin4</th>
      <td>0</td>
      <td>2755325</td>
      <td>1629314</td>
      <td>362552</td>
      <td>16806691</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin1</th>
      <td>0</td>
      <td>1553051</td>
      <td>801500</td>
      <td>197382</td>
      <td>9171735</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin2</th>
      <td>0</td>
      <td>769594</td>
      <td>420527</td>
      <td>94190</td>
      <td>4589512</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin3</th>
      <td>0</td>
      <td>1553359</td>
      <td>870106</td>
      <td>208427</td>
      <td>9511281</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin4</th>
      <td>0</td>
      <td>1719502</td>
      <td>992352</td>
      <td>222066</td>
      <td>10378573</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin1</th>
      <td>0</td>
      <td>1571648</td>
      <td>827637</td>
      <td>202992</td>
      <td>9378229</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin2</th>
      <td>0</td>
      <td>823181</td>
      <td>467869</td>
      <td>108967</td>
      <td>4937604</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin3</th>
      <td>0</td>
      <td>1063439</td>
      <td>531470</td>
      <td>157252</td>
      <td>6489860</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin4</th>
      <td>0</td>
      <td>885233</td>
      <td>483083</td>
      <td>120627</td>
      <td>5328711</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin1</th>
      <td>0</td>
      <td>1582936</td>
      <td>832175</td>
      <td>206263</td>
      <td>9595116</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin2</th>
      <td>0</td>
      <td>685885</td>
      <td>391213</td>
      <td>91151</td>
      <td>4165893</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin3</th>
      <td>0</td>
      <td>1677970</td>
      <td>966261</td>
      <td>224740</td>
      <td>10144858</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin4</th>
      <td>0</td>
      <td>514727</td>
      <td>284784</td>
      <td>66065</td>
      <td>3013631</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin1</th>
      <td>0</td>
      <td>3407768</td>
      <td>1908419</td>
      <td>450279</td>
      <td>20697521</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin2</th>
      <td>0</td>
      <td>1412918</td>
      <td>820282</td>
      <td>187074</td>
      <td>8609649</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin3</th>
      <td>0</td>
      <td>360486</td>
      <td>200904</td>
      <td>46978</td>
      <td>2103701</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin4</th>
      <td>0</td>
      <td>13205</td>
      <td>7707</td>
      <td>1671</td>
      <td>76100</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin1</th>
      <td>0</td>
      <td>4680666</td>
      <td>2672384</td>
      <td>609165</td>
      <td>28335081</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin2</th>
      <td>0</td>
      <td>1027041</td>
      <td>571390</td>
      <td>131322</td>
      <td>6080935</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin3</th>
      <td>0</td>
      <td>99425</td>
      <td>56792</td>
      <td>12921</td>
      <td>578787</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin4</th>
      <td>0</td>
      <td>317</td>
      <td>250</td>
      <td>117</td>
      <td>1968</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin1</th>
      <td>0</td>
      <td>4954380</td>
      <td>2962047</td>
      <td>652872</td>
      <td>29883841</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin2</th>
      <td>0</td>
      <td>194556</td>
      <td>112781</td>
      <td>25183</td>
      <td>1172918</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin3</th>
      <td>0</td>
      <td>3064</td>
      <td>1650</td>
      <td>351</td>
      <td>17234</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin4</th>
      <td>0</td>
      <td>471</td>
      <td>330</td>
      <td>226</td>
      <td>2500</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin1</th>
      <td>0</td>
      <td>5991539</td>
      <td>3580807</td>
      <td>787813</td>
      <td>36170772</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin2</th>
      <td>0</td>
      <td>175983</td>
      <td>104965</td>
      <td>22481</td>
      <td>1081222</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin3</th>
      <td>0</td>
      <td>334</td>
      <td>200</td>
      <td>54</td>
      <td>2096</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin4</th>
      <td>0</td>
      <td>219</td>
      <td>142</td>
      <td>26</td>
      <td>1423</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin1</th>
      <td>0</td>
      <td>5783753</td>
      <td>3342795</td>
      <td>750860</td>
      <td>34862930</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin2</th>
      <td>0</td>
      <td>160069</td>
      <td>87811</td>
      <td>20978</td>
      <td>965723</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin3</th>
      <td>0</td>
      <td>466</td>
      <td>233</td>
      <td>123</td>
      <td>2437</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin4</th>
      <td>0</td>
      <td>286</td>
      <td>146</td>
      <td>12</td>
      <td>1340</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin1</th>
      <td>0</td>
      <td>4914958</td>
      <td>2808349</td>
      <td>639076</td>
      <td>29504829</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin2</th>
      <td>0</td>
      <td>146659</td>
      <td>79970</td>
      <td>18446</td>
      <td>868299</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin3</th>
      <td>0</td>
      <td>371</td>
      <td>224</td>
      <td>90</td>
      <td>2176</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin4</th>
      <td>0</td>
      <td>176</td>
      <td>114</td>
      <td>298</td>
      <td>931</td>
    </tr>
    <tr>
      <th rowspan="40" valign="top">pool2</th>
      <th>SortSeq_bin1</th>
      <td>0</td>
      <td>1535859</td>
      <td>809466</td>
      <td>282203</td>
      <td>6932793</td>
    </tr>
    <tr>
      <th>SortSeq_bin2</th>
      <td>0</td>
      <td>1994122</td>
      <td>1195111</td>
      <td>394137</td>
      <td>9986500</td>
    </tr>
    <tr>
      <th>SortSeq_bin3</th>
      <td>0</td>
      <td>3112859</td>
      <td>2021383</td>
      <td>626649</td>
      <td>16064075</td>
    </tr>
    <tr>
      <th>SortSeq_bin4</th>
      <td>0</td>
      <td>3279820</td>
      <td>2111541</td>
      <td>686484</td>
      <td>17692092</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin1</th>
      <td>0</td>
      <td>1330609</td>
      <td>746110</td>
      <td>269103</td>
      <td>6760551</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin2</th>
      <td>0</td>
      <td>256169</td>
      <td>156675</td>
      <td>53004</td>
      <td>1330057</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin3</th>
      <td>0</td>
      <td>1780228</td>
      <td>1075290</td>
      <td>366477</td>
      <td>9391125</td>
    </tr>
    <tr>
      <th>TiteSeq_01_bin4</th>
      <td>0</td>
      <td>1995933</td>
      <td>1288873</td>
      <td>413947</td>
      <td>10639742</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin1</th>
      <td>0</td>
      <td>1340603</td>
      <td>747697</td>
      <td>272408</td>
      <td>6841566</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin2</th>
      <td>0</td>
      <td>863542</td>
      <td>507085</td>
      <td>176189</td>
      <td>4479884</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin3</th>
      <td>0</td>
      <td>2466965</td>
      <td>1444072</td>
      <td>530148</td>
      <td>12967519</td>
    </tr>
    <tr>
      <th>TiteSeq_02_bin4</th>
      <td>0</td>
      <td>686692</td>
      <td>419611</td>
      <td>144924</td>
      <td>3645739</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin1</th>
      <td>0</td>
      <td>1881279</td>
      <td>1106841</td>
      <td>391059</td>
      <td>9693230</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin2</th>
      <td>0</td>
      <td>1018411</td>
      <td>608125</td>
      <td>205566</td>
      <td>5292042</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin3</th>
      <td>0</td>
      <td>2192254</td>
      <td>1345964</td>
      <td>453446</td>
      <td>11736980</td>
    </tr>
    <tr>
      <th>TiteSeq_03_bin4</th>
      <td>0</td>
      <td>1060345</td>
      <td>657609</td>
      <td>220694</td>
      <td>5587758</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin1</th>
      <td>0</td>
      <td>2515084</td>
      <td>1527876</td>
      <td>522564</td>
      <td>13184739</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin2</th>
      <td>0</td>
      <td>1903439</td>
      <td>1192013</td>
      <td>388293</td>
      <td>10074934</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin3</th>
      <td>0</td>
      <td>1927734</td>
      <td>1215980</td>
      <td>402330</td>
      <td>10254104</td>
    </tr>
    <tr>
      <th>TiteSeq_04_bin4</th>
      <td>0</td>
      <td>158883</td>
      <td>102787</td>
      <td>32322</td>
      <td>809719</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin1</th>
      <td>0</td>
      <td>3836675</td>
      <td>2406011</td>
      <td>796533</td>
      <td>19946374</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin2</th>
      <td>0</td>
      <td>1563516</td>
      <td>969426</td>
      <td>318805</td>
      <td>8221875</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin3</th>
      <td>0</td>
      <td>332584</td>
      <td>205864</td>
      <td>66936</td>
      <td>1700288</td>
    </tr>
    <tr>
      <th>TiteSeq_05_bin4</th>
      <td>0</td>
      <td>11246</td>
      <td>6003</td>
      <td>2335</td>
      <td>58123</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin1</th>
      <td>0</td>
      <td>4534911</td>
      <td>2880466</td>
      <td>934806</td>
      <td>23881666</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin2</th>
      <td>0</td>
      <td>311377</td>
      <td>188856</td>
      <td>62948</td>
      <td>1627496</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin3</th>
      <td>0</td>
      <td>2043</td>
      <td>1321</td>
      <td>449</td>
      <td>11341</td>
    </tr>
    <tr>
      <th>TiteSeq_06_bin4</th>
      <td>0</td>
      <td>3204</td>
      <td>1918</td>
      <td>817</td>
      <td>15389</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin1</th>
      <td>0</td>
      <td>4971451</td>
      <td>3291581</td>
      <td>1035747</td>
      <td>26123287</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin2</th>
      <td>0</td>
      <td>400559</td>
      <td>260676</td>
      <td>81613</td>
      <td>2138543</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin3</th>
      <td>0</td>
      <td>5210</td>
      <td>3167</td>
      <td>1064</td>
      <td>28374</td>
    </tr>
    <tr>
      <th>TiteSeq_07_bin4</th>
      <td>0</td>
      <td>2833</td>
      <td>1795</td>
      <td>532</td>
      <td>15300</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin1</th>
      <td>0</td>
      <td>5631824</td>
      <td>3583466</td>
      <td>1164108</td>
      <td>29642202</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin2</th>
      <td>0</td>
      <td>378038</td>
      <td>233421</td>
      <td>81802</td>
      <td>1994328</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin3</th>
      <td>0</td>
      <td>6144</td>
      <td>3469</td>
      <td>1219</td>
      <td>31033</td>
    </tr>
    <tr>
      <th>TiteSeq_08_bin4</th>
      <td>0</td>
      <td>6557</td>
      <td>3942</td>
      <td>1188</td>
      <td>32173</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin1</th>
      <td>0</td>
      <td>6072547</td>
      <td>3852156</td>
      <td>1262720</td>
      <td>31914552</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin2</th>
      <td>0</td>
      <td>247557</td>
      <td>151018</td>
      <td>51124</td>
      <td>1300217</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin3</th>
      <td>0</td>
      <td>287</td>
      <td>212</td>
      <td>90</td>
      <td>1373</td>
    </tr>
    <tr>
      <th>TiteSeq_09_bin4</th>
      <td>0</td>
      <td>92</td>
      <td>77</td>
      <td>334</td>
      <td>405</td>
    </tr>
  </tbody>
</table>


Now we plot the barcode-read fates for each library / sample, showing the bars for valid barcodes in orange and the others in gray.
We see that the largest fraction of barcode reads correspond to valid barcodes, and most of the others are invalid barcodes (probably because the map to variants that aren't present in our variant table since we didn't associate all variants with barcodes). The exception to this is lib2 Titeseq_03_bin3; the PCR for this sample in the original sequencing run failed, so we followed it up with a single MiSeq lane. We did not filter out the PhiX reads from this data before parsing, so these PhiX reads will deflate the fraction of valid barcode reads as expected, but does not indicate any problems.


```python
barcode_fate_plot = (
    ggplot(
        fates
        .assign(sample=lambda x: pd.Categorical(x['sample'],
                                                x['sample'].unique(),
                                                ordered=True),
                fate=lambda x: pd.Categorical(x['fate'],
                                              x['fate'].unique(),
                                              ordered=True),
                is_valid=lambda x: x['fate'] == 'valid barcode'
                ), 
        aes('fate', 'count', fill='is_valid')) +
    geom_bar(stat='identity') +
    facet_grid('sample ~ library') +
    facet_grid('sample ~ library') +
    scale_fill_manual(CBPALETTE, guide=False) +
    theme(figure_size=(1.4 * (1 + fates['library'].nunique()),
                       1.7 * (1.2 + fates['sample'].nunique())),
          axis_text_x=element_text(angle=90),
          panel_grid_major_x=element_blank()
          ) +
    scale_y_continuous(labels=dms_variants.utils.latex_sci_not,
                       name='number of reads')
    )

_ = barcode_fate_plot.draw()
```


    
![png](count_variants_files/count_variants_42_0.png)
    


## Output csv of barcode counts in variant-barcode lookup table


```python
print(f"Writing variant counts to {config['variant_counts_file']}")
counts.to_csv(config['variant_counts_file'], index=False)
```

    Writing variant counts to results/counts/variant_counts.csv


The [CodonVariantTable](https://jbloomlab.github.io/dms_variants/dms_variants.codonvarianttable.html#dms_variants.codonvarianttable.CodonVariantTable) has lots of nice functions that can be used to analyze the counts it contains.
However, we do that in the next notebook so we don't have to re-run this entire (rather computationally intensive) notebook every time we want to analyze a new aspect of the counts.


```python

```
