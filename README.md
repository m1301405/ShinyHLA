# ShinyHLA: an R/Shiny-based web application for HLA typing of NGS data
ShinyHLA is an R/Shiny-based tool for HLA typing, supporting WES and RNA-seq data, integrating leading tools and databases, and offering a user-friendly interface with advanced visualization and analysis features.

- *Data compatibility: Supports HLA typing for both WES and RNA-seqdata.*
- *Integration of tools and database: Incorporates OptiType, arcasHLA, HLA-HD, and SpecHLA, along with the latest IPD-IMGT/HLA database.*
- *User-friendly interface: Provides a graphical user interface (GUI) that eliminates the need for programming experience.*
- *Enhanced data visualization and analysis: Incorporates a pivot table and the Integrative Genomics Viewer (IGV) for efficient result review and comparison.*

## Website
The official web documentation for ShinyHLA are available online. 
- Documentation: [ShinyHLA documentation](https://animated-reminder-9ff.notion.site/ebd//7b6f9fcf9d034442b2b64c89f8662a6c)

It provides:
1.	Step-by-step guides: Detailed instructions for tasks such as file upload, preprocessing, quality control, and HLA typing.
2.	Tutorials: Clear walkthroughs for both beginners and advanced users, demonstrating core features like HLA typing, result visualization using IGV, and pivot table usage for data integration.
3.	Docker deployment instructions: Simple and clear guidelines for setting up Docker environments to optimize preprocessing and streamline the analysis workflow.

## Requirements
ShinyHLA uses the following software and reference: 
### 1. HLA typing tools
- [OptiType v1.3.5](https://github.com/FRED-2/OptiType)
- [arcasHLA v0.6.0](https://github.com/RabadanLab/arcasHLA)
- [HLA-HD v1.7.0](https://w3.genome.med.kyoto-u.ac.jp/HLA-HD/)
- [SpecHLA v1.0.7](https://github.com/deepomicslab/SpecHLA)

### 2. R packages
- Please refer to [ShinyHLA documentation](https://animated-reminder-9ff.notion.site/ebd//7b6f9fcf9d034442b2b64c89f8662a6c), section 1.7.

### 3. Sequence Processing Tools
- [seqtk v1.3-r106](https://github.com/lh3/seqtk)
- [Pandepth v2.25](https://github.com/HuiyangYu/PanDepth)

### 4. Human genome reference
- [Gencode](https://www.gencodegenes.org/human/release_38.html): gencode.v38.annotation.gff3

## Installation via Docker 
1. Install Docker on your computer and make sure it works.
2. Call `docker pull xxxx/shinyhla` which will download the Docker image.
3. You can use the image as followes: \
   `docker run -d -p 3838:3838 --name shinyhla shinyhla:latest`
4. Open your web browser and navigate to `http://localhost:3838`.

For more detailed installation and operation instructions, please refer to [ShinyHLA documentation](https://animated-reminder-9ff.notion.site/ebd//7b6f9fcf9d034442b2b64c89f8662a6c), section 3.

## Citation
[Cheng-Chi Yang, xxx, xxx, and xxx. "ShinyHLA: An R/Shiny-Based Web Application for HLA Typing of NGS Data." journal (Date).](https://unknown.com)

## Getting help
Should you have any queries, please feel free to contact us, we will reply as soon as possible ([m1301405@cgu.edu.tw](mailto:m1301405@cgu.edu.tw)).



