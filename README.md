# UManSysProp_public
This repository presents the source code for predictive techniques available in the UManSysProp facility.  Please note these files alone will not enable you to run your own version of the UManSysProp server as demonstrated by our website. Rather, the files provided here are to enable you to call our functions from your own Python modules. To run your own server, please download the approproate package from https://github.com/waveform80/umansysprop and then insert the files and folders provided here within the umansysprop sub-directory.

All files presented here are covered under the GNU GPL license v3.0. For more information, please read the license file. A brief overview of the viable permissions can be found here: http://choosealicense.com/licenses/

The UManSysProp facility was created to provide automated predictions of molecular and mixture properties. The website  http://umansysprop.seaes.manchester.ac.uk/ provides details on the provenance of the predictive techniques provided here, as well as providing you with a web based interface for carrying out predictions without handling any code. For more flexibility, we also provide a programmer friendly JSON API that enables you to call our suite of tools from your own code without opening a web browser. More information can be found on our site.  

To re-iterate, here we provide the source code for those predictive techniques outside of the server framework. Written in Python, the UManSysProp folder can be imported as a module which then gives relevant access to each technique. For example, to calculate the functional groups for use with the AIOMFAC technique, one might write:
<br />
import umansysprop.groups #umansysprop is constructed as a package <br />
import pybel <br />
AIOMFAC_keys = {}<br />
for s in SMILES:<br />
    SMILES_object=pybel.readstring(b'smi',s)<br />
    AIOMFAC_keys[s] = umansysprop.groups.aiomfac(SMILES_object)<br />

Note you need to have the Python interface to the openbabel package installed (Pybel).
<br />
Ultimately, we aim to generate a user community around UManSysProp. As we develop and evaluate more predictive techniques and models, according to the scientific peer review process, these will be included in future releases. If you use any of our developments in your code, following the license restrictions, we would appreciate hearing about it. We also request you reference our development paper:
<br />
[Topping, D., Barley, M., Bane, M. K., Higham, N., Aumont, B., Dingle, N., and McFiggans, G.: UManSysProp v1.0: an online and open-source facility for molecular property prediction and atmospheric aerosol calculations, Geosci. Model Dev., 9, 899-914, doi:10.5194/gmd-9-899-2016, 2016.](http://www.geosci-model-dev.net/9/899/2016/)
<br />
Code DOI [![DOI](https://zenodo.org/badge/20123/loftytopping/UManSysProp_public.svg)](https://zenodo.org/badge/latestdoi/20123/loftytopping/UManSysProp_public)



In addition, if you would like any particular technique added to this facility then get in touch and we can discuss.

For up to date information, including paper and training alerts, please check our website:http://umansysprop.seaes.manchester.ac.uk/ . You can also follow us on twitter @UManSysProp

