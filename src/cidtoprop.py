# quickprops/modules/cid_to_properties.py
import requests

def fetch_pubchem_properties(cid):
    """
    Fetch properties from PubChem using CID.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/MolecularFormula,MolecularWeight,IUPACName,CanonicalSMILES/JSON"
    response = requests.get(url)
    props = response.json()['PropertyTable']['Properties'][0]
    return props

