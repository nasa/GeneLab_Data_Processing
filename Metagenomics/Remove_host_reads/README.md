# GeneLab pipeline for removing host reads in metagenomics sequencing data

> **It is NASA's policy that any human reads are to be removed from metagenomics datasets prior to being hosted in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/). As such, the document [`GL-DPPD-7105-B.md`](Pipeline_GL-DPPD-7105_Versions/GL-DPPD-7105-B.md) holds an overview and example commands for how GeneLab identifies and removes human DNA in metagenomics sequencing datasets. See the [Repository Links](#repository-links) descriptions below for more information. The percentage of human reads removed and a GeneLab human read removal summary is provided for each GLDS dataset in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).**  
>
>**Because host-read removal is broadly needed for metagenomics processing, the original _MGRemoveHumanReads_ has been expanded to include removal of other host DNA and renamed _MGRemoveHostReads_. The current pipeline supports removal of human reads by default as well as reads from any other host organism relevant to the dataset.** 
> 
> Note: The exact human read identification and removal commands as well as pipeline version used for specific GLDS datasets can be found in the *_processing_info.zip file under "Files" for each respective GLDS dataset in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/). 
---

## Repository Links

* [**Pipeline_GL-DPPD-7105_Versions**](Pipeline_GL-DPPD-7105_Versions)

  - Contains the versions documentation of both the current GeneLab pipeline for identifying and removing host reads in metagenomics sequencing data (MGRemoveHostReads) and the previous GeneLab pipeline dedicated to identifying and removing human reads only (MGRemoveHumanReads)

* [**Workflow_Documentation**](Workflow_Documentation)

  - Contains instructions for installing and running the current GeneLab MGRemoveHostReads workflow and the previous GeneLab MGRemoveHumanReads workflow

---

**Developed and maintained by:**  
Michael D. Lee (Mike.Lee@nasa.gov)  
Jihan Yehia (jihan.yehia@nasa.gov)
