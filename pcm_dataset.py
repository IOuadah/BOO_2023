from qsprpred.data.sources.papyrus import Papyrus
from qsprpred.models.tasks import TargetTasks
from qsprpred.data.data import QSPRDataset
from qsprpred.extra.data.data import PCMDataset



def AR_PCM_both(activity, dataset_name, data_dir='data'):
    """
    A classification dataset that contains activity data for a PCM approach to model activity for a selection of adenosine receptors. The function recreates steps from data_preparation_advanced.ipynb.

    Returns:
        a `QSPRDataset` instance with the loaded data
    """

    acc_keys = ["P29274", "P29275", "P30542", "P0DMS8", "P25099", "P30543", "P29276", "P28647"]
    quality = "high"  # choose minimum quality from {"high", "medium", "low"}
    papyrus_version = '05.6'  # Papyrus database version

    papyrus = Papyrus(
        data_dir=data_dir,
        stereo=False,
        version=papyrus_version
    )

    mt = papyrus.getData(
        acc_keys,
        quality,
        name=dataset_name,
        use_existing=True,
        activity_types=activity
    )
    
    mt.standardizeSmiles('chembl', drop_invalid=True)
    # mt.standardizeSmiles('chembl', drop_invalid=False)
    mt.dropInvalids()

    ds_seq = papyrus.getProteinData(acc_keys, name=f"{mt.name}_seqs", use_existing=True)

    def sequence_provider(acc_keys):
        """
        A function that provides a mapping from accession key to a protein sequence.

        Args:
            acc_keys (list): Accession keys of the protein to get a mapping of sequences for.

        Returns:
            (dict) : Mapping of accession keys to protein sequences.
            (dict) : Additional information to pass to the MSA provider (can be empty).
        """
        map = dict()
        info = dict()
        for i, row in ds_seq.iterrows():
            map[row['accession']] = row['Sequence']

            # can be omitted
            info[row['accession']] = {
                'Organism': row['Organism'],
                'UniProtID': row['UniProtID'],
            }

        return map, info

    return PCMDataset.fromMolTable(
        mt,
        target_props=[
            {
                "name": "pchembl_value_Median",
                "task": TargetTasks.SINGLECLASS,
                "th": [6.5]
            }
        ],
        proteincol="accession",
        proteinseqprovider=sequence_provider,
    )

def AR_PCM_HUMAN(activity, dataset_name, data_dir='data'):
    """
    A classification dataset that contains activity data for a PCM approach to model activity for a selection of adenosine receptors. The function recreates steps from data_preparation_advanced.ipynb.

    Returns:
        a `QSPRDataset` instance with the loaded data
    """

    acc_keys = ["P29274", "P29275", "P30542", "P0DMS8"]
    quality = "high"  # choose minimum quality from {"high", "medium", "low"}
    papyrus_version = '05.6'  # Papyrus database version

    papyrus = Papyrus(
        data_dir=data_dir,
        stereo=False,
        version=papyrus_version
    )

    mt = papyrus.getData(
        acc_keys,
        quality,
        name=dataset_name,
        use_existing=True,
        activity_types=activity
    )
    
    mt.standardizeSmiles('chembl', drop_invalid=True)
    # mt.standardizeSmiles('chembl', drop_invalid=False)
    mt.dropInvalids()
    
    ds_seq = papyrus.getProteinData(acc_keys, name=f"{mt.name}_seqs", use_existing=True)

    def sequence_provider(acc_keys):
        """
        A function that provides a mapping from accession key to a protein sequence.

        Args:
            acc_keys (list): Accession keys of the protein to get a mapping of sequences for.

        Returns:
            (dict) : Mapping of accession keys to protein sequences.
            (dict) : Additional information to pass to the MSA provider (can be empty).
        """
        map = dict()
        info = dict()
        for i, row in ds_seq.iterrows():
            map[row['accession']] = row['Sequence']

            # can be omitted
            info[row['accession']] = {
                'Organism': row['Organism'],
                'UniProtID': row['UniProtID'],
            }

        return map, info

    return PCMDataset.fromMolTable(
        mt,
        target_props=[
            {
                "name": "pchembl_value_Median",
                "task": TargetTasks.SINGLECLASS,
                "th": [6.5]
            }
        ],
        proteincol="accession",
        proteinseqprovider=sequence_provider,
    )

def AR_PCM_RAT(activity, dataset_name, data_dir='data'):
    """
    A classification dataset that contains activity data for a PCM approach to model activity for a selection of adenosine receptors. The function recreates steps from data_preparation_advanced.ipynb.

    Returns:
        a `QSPRDataset` instance with the loaded data
    """

    acc_keys = ["P25099", "P30543", "P29276", "P28647"]
    quality = "high"  # choose minimum quality from {"high", "medium", "low"}
    papyrus_version = '05.6'  # Papyrus database version

    papyrus = Papyrus(
        data_dir=data_dir,
        stereo=False,
        version=papyrus_version
    )

    mt = papyrus.getData(
        acc_keys,
        quality,
        name=dataset_name,
        use_existing=False,
        activity_types=activity
    )
    
    mt.standardizeSmiles('chembl', drop_invalid=True)
    # mt.standardizeSmiles('chembl', drop_invalid=False)
    mt.dropInvalids()
    
    ds_seq = papyrus.getProteinData(acc_keys, name=f"{mt.name}_seqs", use_existing=True)

    def sequence_provider(acc_keys):
        """
        A function that provides a mapping from accession key to a protein sequence.

        Args:
            acc_keys (list): Accession keys of the protein to get a mapping of sequences for.

        Returns:
            (dict) : Mapping of accession keys to protein sequences.
            (dict) : Additional information to pass to the MSA provider (can be empty).
        """
        map = dict()
        info = dict()
        for i, row in ds_seq.iterrows():
            map[row['accession']] = row['Sequence']

            # can be omitted
            info[row['accession']] = {
                'Organism': row['Organism'],
                'UniProtID': row['UniProtID'],
            }

        return map, info

    return PCMDataset.fromMolTable(
        mt,
        target_props=[
            {
                "name": "pchembl_value_Median",
                "task": TargetTasks.SINGLECLASS,
                "th": [6.5]
            }
        ],
        proteincol="accession",
        proteinseqprovider=sequence_provider,
    )