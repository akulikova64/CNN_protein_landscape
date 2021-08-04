import argparse
from pathlib import Path
from wget import download, bar_thermometer


def cli():

    parser = argparse.ArgumentParser()

    parser.add_argument('pdb_file', type=Path)

    parser.add_argument('-d','--data-dir', required=False, default=Path("./pdbs"), type=Path)

    parser.add_argument('-t', '--pdb-type', required=False, choices=['pdb','cif', 'both'], default='both')

    parser.add_argument('--api', required=False, choices=['rcsb','pdb-redo'], default='rcsb')

    args = parser.parse_args()

    args.data_dir.mkdir(0o744, parents=True, exist_ok=True)

    return args



def download_proteins(pdb_file: Path, data_dir: Path, pdb_type='both', api='rcsb'):

    import pymp

    url = '' 
    if api == 'rcsb':
        url = 'http://files.rcsb.org/download/{}'

    elif api == 'pdb-redo':
        url = 'https://pdb-redo.eu/db/{}/{}'

    if pdb_type is 'both':
        pdb_type = ['pdb', 'cif']

    with pdb_file.open('rt') as f:
        pdb_list = [line.strip() for line in f.readlines() if not '#' in line]

    total_proteins = len(pdb_list)

    failed_proteins = pymp.shared.list()
    with pymp.Parallel() as p:
        for idx in p.xrange(total_proteins):
            pdb_id = pdb_list[idx]
            try: 
                # TODO refactor with requests?
                if api == 'rcsb':
                    if 'pdb' in pdb_type:
                        download(
                            url.format(pdb_id.upper() + '.pdb'), 
                            out=str(data_dir/ (pdb_id.lower() + '.pdb')),
                            bar=bar_thermometer
                        )

                    if 'cif' in pdb_type:
                        download(
                            url.format(pdb_id.upper() + '.cif'), 
                            out=str(data_dir/ (pdb_id.lower() + '.cif')),
                            bar=bar_thermometer    
                        )

                elif api == 'pdb-redo':
                    if 'pdb' in pdb_type:
                        filename = pdb_id + '_final_tot.pdb'
                        download(
                            url.format(pdb_id, filename), 
                            out=str(data_dir / filename),
                            bar=bar_thermometer    
                        )

                    if 'cif' in pdb_type:
                        filename = pdb_id + '_final.cif'
                        download(
                            url.format(pdb_id, filename), 
                            out=str(data_dir / filename),
                            bar=bar_thermometer    
                        )

            except Exception as e:
                failed_proteins.append(pdb_id)
                p.print(f'Error with protein ({pdb_id}):\n{e}\n')

            finally:
                p.print(f'Thread: {p.threadcd um} Finished protein ({pdb_id}): {idx+1} out of {total_proteins}')

    log_failed_proteins = data_dir / 'failed_pdbs.txt'

    if len(failed_proteins) > 0:
        with log_failed_proteins.open('wt') as f:
            f.write(f'#Dataset: {pdb_file}\n')
            f.write('\n'.join(failed_proteins))



if __name__=="__main__":
    '''$ PYMP_NUM_THREADS=<number_of_threads> python download_proteins.py path/to/pdb_dataset.txt'''
    args = cli()

    download_proteins(**vars(args))
