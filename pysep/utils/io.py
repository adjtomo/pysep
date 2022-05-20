"""
Read utilities for Pysep
"""
import re
import yaml
import numpy as np


def read_yaml(fid):
    """

    :return:
    """
    # work around PyYAML bugs
    yaml.SafeLoader.add_implicit_resolver(
        u'tag:yaml.org,2002:float',
        re.compile(u'''^(?:
         [-+]?(?:[0-9][0-9_]*)\\.[0-9_]*(?:[eE][-+]?[0-9]+)?
        |[-+]?(?:[0-9][0-9_]*)(?:[eE][-+]?[0-9]+)
        |\\.[0-9_]+(?:[eE][-+][0-9]+)?
        |[-+]?[0-9][0-9_]*(?::[0-5]?[0-9])+\\.[0-9_]*
        |[-+]?\\.(?:inf|Inf|INF)
        |\\.(?:nan|NaN|NAN))$''', re.X),
        list(u'-+0123456789.'))

    with open(fid, 'r') as f:
        config = yaml.safe_load(f)

    # Replace 'None' and 'inf' values to match expectations
    for key, val in config.items():
        if val == "None":
            config[key] = None
        if val == "inf":
            config[key] = np.inf

    return config
