# quickprops/modules/smiles_to_cid.py
import requests

def get_pubchem_cid(smiles):
    """
    Get PubChem CID from SMILES.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/cids/JSON"
    response = requests.get(url)
    data = response.json()
    cids = data.get('IdentifierList', {}).get('CID', [])
    if not cids:
        return None
    return cids[0]

