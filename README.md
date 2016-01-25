# UManSysProp_public
This repository presents the source code for predictive techniques available in the UManSysProp facility.  All files presented here are covered under the GNU GPL license v3.0. For more information, please read the license file. A brief overview of the viable permissions can be found here: http://choosealicense.com/licenses/

The UManSysProp facility was created to provide automated predictions of molecular and mixture properties. The website  http://umansysprop.seaes.manchester.ac.uk/ provides details on the provenance of the predictive techniques provided, as well as providing you with a web based interface for carrying out predictions without handling any code. For more flexibility, we also provide a programmer friendly JSON API that enables you to call our suite of tools from your own code without opening a web browser. More information can be found on our site.  Here we provide the source code for those predictive techniques. Written in Python, the UManSysProp folder can be imported as a module which then gives relevant access to each technique. For example, to calculate the functional groups for use with the AIOMFAC technique, one might write:
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
[Topping, D., Barley, M. H., Bane, M., Higham, N., Aumont, B., and McFiggans, G.: UManSysProp: an online facility for molecular property prediction and atmospheric aerosol calculations, Geosci. Model Dev. Discuss., 8, 9669-9706, doi:10.5194/gmdd-8-9669-2015, 2015.](http://www.geosci-model-dev-discuss.net/gmd-2015-197/)
<br />
In addition, if you would like any particular technique added to this facility then get in touch and we can discuss.

For up to date information, including paper and training alerts, please check our website:http://umansysprop.seaes.manchester.ac.uk/ . You can also follow us on twitter @UManSysProp

