# ----------------------------------------------------------------------------
# Copyright (c) 2011-2015, The BIOM Format Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import sys
import biom
import numpy as np
import tempfile
import h5py

if len(sys.argv) < 3:
    raise SystemExit

if '://' in sys.argv[1]:
    from urllib import request
    fp, _ = request.urlretrieve(sys.argv[1])
else:
    fp = sys.argv[1]
axis = sys.argv[2]
table = biom.load_table(fp)

# drop some samples, verify what gets dropped
sample_sums = table.sum(axis=axis)
idx = np.argsort(sample_sums)[:int(len(table.ids()) / 2)]
dropped_ids = table.ids(axis=axis)[idx]
obs_by_id = table.filter(lambda v, i, md: i not in set(dropped_ids),
                         axis=axis, inplace=False)
assert set(obs_by_id.ids()).intersection(set(dropped_ids)) == set()
obs_by_id = table.filter(set(dropped_ids), invert=True,
                         axis=axis, inplace=False)
assert set(obs_by_id.ids()).intersection(set(dropped_ids)) == set()

obs_without_top = table.filter(lambda v, i, md: v.sum() < sample_sums.max(),
                               inplace=False, axis=axis)
exp = table.ids(axis=axis)[sample_sums != sample_sums.max()]
assert set(obs_without_top.ids(axis=axis)) == set(exp)

# arbitrary partitioning and collapsing
md = {i: {'partition': True if idx % 2 else False}
      for idx, i in enumerate(table.ids(axis=axis))}
table.add_metadata(md, axis=axis)
parts = list(table.partition(lambda i, m: m['partition'], axis=axis))
collapsed = table.collapse(lambda i, m: m['partition'], axis=axis, norm=False)

assert len(parts) == 2
assert (len(parts[0][1].ids(axis=axis)) + len(parts[1][1].ids(axis=axis))) == \
        len(table.ids(axis=axis))
collapsed_sums = collapsed.sum(axis=axis)
parts_as_dict = dict(parts)
for name, obs in zip(collapsed.ids(axis=axis), collapsed_sums):
    assert parts_as_dict[name].sum() == obs

regrouped = parts[0][1].concat([parts[1][1]], axis=axis)
regrouped = regrouped.sort_order(table.ids(axis=axis), axis=axis)
invaxis = 'sample' if axis == 'observation' else 'observation'
regrouped = regrouped.sort_order(table.ids(axis=invaxis), axis=invaxis)
assert regrouped == table

regrouped.del_metadata(keys=['partition'], axis=axis)
assert regrouped != table
table.del_metadata(keys=['partition'], axis=axis)
assert regrouped == table

# transforms
pa = table.pa(inplace=False)
assert pa.sum() == table.matrix_data.nnz


class inc:
    def __init__(self):
        self.scaler = 0

    def __call__(self, v, i, md):
        new_v = v + self.scaler
        self.scaler += 1
        return new_v


scaled = pa.transform(inc(), inplace=False, axis=axis)
for idx, v in enumerate(scaled.iter_data(axis=axis, dense=False)):
    assert np.allclose(v.data, np.ones(1) + idx)

# roundtrips
_, filepath = tempfile.mkstemp()
with h5py.File(filepath, 'w') as out:
    table.to_hdf5(out, 'rc-testing')
obs = biom.load_table(filepath)
assert obs == table

with open(filepath, 'w') as out:
    out.write(table.to_json('rc-testing'))
obs = biom.load_table(filepath)
assert obs == table

with open(filepath, 'w') as out:
    out.write(table.to_tsv())
obs = biom.load_table(filepath)
exp = table.copy()

# TSV metadata are a nightmare especially as there are custom parsing
# is necessary for things like taxonomy.
exp.del_metadata()
exp.type = None  # a TSV does not have a table type
assert obs == exp
