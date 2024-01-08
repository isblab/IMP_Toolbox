import IMP
import RMF
import IMP.atom
import IMP.core
import IMP.rmf
import IMP.pmi.analysis
import sys
import re
import tqdm


def extract_names_and_sort(particles_rmf, particles_pdb):
    pdb_names = []
    for n in particles_pdb:
        name = n.get_name()
        match = re.search('Atom CA of residue [0-9]{3}', name)
        assert match, f'Regex failed for pdb particle: {name}'
        resnum = int(name[-3:])
        pdb_names.append(resnum)
    rmf_names = []
    rmf_selected = []
    for i, n in enumerate(particles_rmf):
        name = n.get_name()
        if '_bead' in name:
            assert not ('Fragment' in name)
            rmf_names.append(name)
            continue  # skip the flexible beads (not in pdb)
        match = re.search('[0-9]{3}-[0-9]{3}', name)
        assert match, f'Regex failed for rmf particle: {name}'
        resnum = (int(match.group(0)[:3]), int(match.group(0)[-3:]))
        resnum = resnum[0] + (resnum[1] - resnum[0]) // 2
        rmf_names.append(resnum)
        if resnum in pdb_names:
            rmf_selected.append(i)

    pdb_selected = [pdb_names.index(rmf_names[i]) for i in rmf_selected]
    return rmf_names, rmf_selected, pdb_names, pdb_selected
    

def find_transform_and_save(query, template, output_name, pdb_file, chain_selector):
    if not isinstance(query, dict):
        assert not isinstance(template, dict)
        query = {'A': query}
        template = {'A': template}
        
    new_mdl = IMP.Model()
    reload = IMP.atom.read_pdb(pdb_file, new_mdl, chain_selector)
    tra, trb = IMP.pmi.analysis.Alignment(query=query, template=template).align()
    IMP.atom.transform(reload, trb)
    IMP.atom.write_pdb(reload, output_name)


def extract_coordinates(ccm_file, pdb_file, chain_id, protein_name, copy_index):
    ccm_mdl = IMP.Model()
    ccm = RMF.open_rmf_file_read_only(ccm_file)
    hier = IMP.rmf.create_hierarchies(ccm, ccm_mdl)[0]
    IMP.rmf.load_frame(ccm, 0)
    ccm_mdl.update()
    
    pdb_ca_mdl = IMP.Model()
    pdb_ca = IMP.atom.read_pdb(pdb_file, pdb_ca_mdl, IMP.atom.CAlphaPDBSelector())
    pdb_ca_mdl.update()
    
    sel_pdb_ca = IMP.atom.Selection(pdb_ca, chain_id=chain_id).get_selected_particles()
    sel_ccm = IMP.atom.Selection(hier, molecule=protein_name, copy_index=copy_index).get_selected_particles()
    assert all([IMP.atom.Fragment.get_is_setup(x) for x in sel_ccm]), 'Something is wrong'
    
    rmf_names, rmf_selected, pdb_names, pdb_selected = extract_names_and_sort(sel_ccm, sel_pdb_ca)

    coords_rmf = [IMP.core.XYZ(sel_ccm[i]).get_coordinates() for i in rmf_selected]
    coords_pdb = [IMP.core.XYZ(sel_pdb_ca[i]).get_coordinates() for i in pdb_selected]
    return coords_rmf, coords_pdb
    

ccm_file = sys.argv[1]




# align PKPs
print('Aligning PKP1a')
pdb_file = '1xm9.pdb'
for i in tqdm.trange(4):
    template, query = extract_coordinates(ccm_file, pdb_file, 'A', 'PKP1a', i)
    find_transform_and_save(query, template, f'aligned_pkp1a_{i}.pdb', pdb_file, IMP.atom.ChainPDBSelector('A'))

# align DPs
print('Aligning DP')
pdb_file = '3r6n.pdb'
for i in tqdm.trange(4):
    template, query = extract_coordinates(ccm_file, pdb_file, 'A', 'DP', i)
    find_transform_and_save(query, template, f'aligned_dp_{i}.pdb', pdb_file, IMP.atom.ChainPDBSelector('A'))


# align PG-DSCs
print('Aligning PG-DSCs')
pdb_file = 'pg-dsc1-pg-dsc1_modified.pdb'
for i in tqdm.trange(2):
    template1, query1 = extract_coordinates(ccm_file, pdb_file, 'A', 'PG', i)
    template2, query2 = extract_coordinates(ccm_file, pdb_file, 'C', 'DSC1', i)
    find_transform_and_save({'A': query1, 'B': query2}, {'A': template1, 'B': template2}, f'aligned_pg_dsc1_{i}.pdb', pdb_file,
                            IMP.atom.ChainPDBSelector('A') | IMP.atom.ChainPDBSelector('C'))

# align PG-DSGs
print('Aligning PG-DSGs')
pdb_file = 'pg-dsg1-pg-dsg1_modified.pdb'
for i in tqdm.trange(2):
    template1, query1 = extract_coordinates(ccm_file, pdb_file, 'A', 'PG', i + 2)
    template2, query2 = extract_coordinates(ccm_file, pdb_file, 'C', 'DSG1', i)
    find_transform_and_save({'A': query1, 'B': query2}, {'A': template1, 'B': template2}, f'aligned_pg_dsg1_{i}.pdb', pdb_file,
                            IMP.atom.ChainPDBSelector('A') | IMP.atom.ChainPDBSelector('C'))