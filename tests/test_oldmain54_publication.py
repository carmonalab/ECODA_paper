"""Finite fixed54 metadata/copy boundary; all matrices are tiny synthetic fixtures."""
import copy,hashlib,importlib.util,json,os,pathlib,subprocess,sys,tempfile,unittest
from types import SimpleNamespace
ROOT=pathlib.Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('publisher54',ROOT/'src/utils/py/bl_science_publication.py');pub=importlib.util.module_from_spec(spec);spec.loader.exec_module(pub)
# Public method/parameter/filename scope only; no patient identifiers or scientific payload.
FIXTURE={'batch_effect/corrected_final': {'core': ['embeddings/Alzheimer_batch_effect_corrected_final_hvg2000_highres_mrvi_dists.feather',
                                           'embeddings/Alzheimer_batch_effect_corrected_final_hvg2000_highres_pilot_dists.feather',
                                           'embeddings/Breast_cancer_batch_effect_corrected_final_hvg2000_highres_mrvi_dists.feather',
                                           'embeddings/Breast_cancer_batch_effect_corrected_final_hvg2000_highres_pilot_dists.feather',
                                           'embeddings/Covid19_PBMC_batch_effect_corrected_final_hvg2000_highres_mrvi_dists.feather',
                                           'embeddings/Kidney_KPMP_full_batch_effect_corrected_final_hvg2000_highres_mrvi_dists.feather',
                                           'embeddings/Stephenson_batch_effect_corrected_final_hvg2000_highres_mrvi_dists.feather',
                                           'embeddings/Stephenson_batch_effect_corrected_final_hvg2000_highres_pilot_dists.feather'],
                                  'rows': [{'dataset': 'Stephenson',
                                            'method': 'mrvi',
                                            'parameter': '',
                                            'view': 'batch_effect_corrected'},
                                           {'dataset': 'Alzheimer',
                                            'method': 'mrvi',
                                            'parameter': '',
                                            'view': 'batch_effect_corrected'},
                                           {'dataset': 'Breast_cancer',
                                            'method': 'mrvi',
                                            'parameter': '',
                                            'view': 'batch_effect_corrected'},
                                           {'dataset': 'Covid19_PBMC',
                                            'method': 'mrvi',
                                            'parameter': '',
                                            'view': 'batch_effect_corrected'},
                                           {'dataset': 'Kidney_KPMP_full',
                                            'method': 'mrvi',
                                            'parameter': '',
                                            'view': 'batch_effect_corrected'},
                                           {'dataset': 'Stephenson',
                                            'method': 'pilot',
                                            'parameter': '',
                                            'view': 'batch_effect_corrected'},
                                           {'dataset': 'Alzheimer',
                                            'method': 'pilot',
                                            'parameter': '',
                                            'view': 'batch_effect_corrected'},
                                           {'dataset': 'Breast_cancer',
                                            'method': 'pilot',
                                            'parameter': '',
                                            'view': 'batch_effect_corrected'}],
                                  'runtime': ['embeddings/Alzheimer_batch_effect_corrected_final_hvg2000_highres_mrvi_dists.feather.runtime.json',
                                              'embeddings/Alzheimer_batch_effect_corrected_final_hvg2000_highres_pilot_dists.feather.runtime.json',
                                              'embeddings/Breast_cancer_batch_effect_corrected_final_hvg2000_highres_mrvi_dists.feather.runtime.json',
                                              'embeddings/Breast_cancer_batch_effect_corrected_final_hvg2000_highres_pilot_dists.feather.runtime.json',
                                              'embeddings/Covid19_PBMC_batch_effect_corrected_final_hvg2000_highres_mrvi_dists.feather.runtime.json',
                                              'embeddings/Kidney_KPMP_full_batch_effect_corrected_final_hvg2000_highres_mrvi_dists.feather.runtime.json',
                                              'embeddings/Stephenson_batch_effect_corrected_final_hvg2000_highres_mrvi_dists.feather.runtime.json',
                                              'embeddings/Stephenson_batch_effect_corrected_final_hvg2000_highres_pilot_dists.feather.runtime.json'],
                                  'timings': ['execution_times_batch_effect_corrected_final_mrvi_Alzheimer.feather',
                                              'execution_times_batch_effect_corrected_final_mrvi_Breast_cancer.feather',
                                              'execution_times_batch_effect_corrected_final_mrvi_Covid19_PBMC.feather',
                                              'execution_times_batch_effect_corrected_final_mrvi_Kidney_KPMP_full.feather',
                                              'execution_times_batch_effect_corrected_final_mrvi_Stephenson.feather',
                                              'execution_times_batch_effect_corrected_final_pilot_Alzheimer.feather',
                                              'execution_times_batch_effect_corrected_final_pilot_Breast_cancer.feather',
                                              'execution_times_batch_effect_corrected_final_pilot_Stephenson.feather']},
 'batch_effect/uncorrected_final': {'core': ['embeddings/Alzheimer_batch_effect_uncorrected_final_hvg2000_highres_mrvi_dists.feather',
                                             'embeddings/Breast_cancer_batch_effect_uncorrected_final_hvg2000_highres_mrvi_dists.feather',
                                             'embeddings/Joanito_batch_effect_uncorrected_final_hvg2000_highres_mrvi_dists.feather',
                                             'embeddings/Stephenson_batch_effect_uncorrected_final_hvg2000_highres_mrvi_dists.feather',
                                             'results/Alzheimer_batch_effect_uncorrected_final_gloscope.rds',
                                             'results/Breast_cancer_batch_effect_uncorrected_final_gloscope.rds',
                                             'results/Covid19_PBMC_batch_effect_uncorrected_final_gloscope.rds',
                                             'results/Joanito_batch_effect_uncorrected_final_gloscope.rds',
                                             'results/Kidney_KPMP_full_batch_effect_uncorrected_final_gloscope.rds',
                                             'results/Lung_batch_effect_uncorrected_final_gloscope.rds',
                                             'results/Lupus_PBMC_batch_effect_uncorrected_final_gloscope.rds',
                                             'results/Stephenson_batch_effect_uncorrected_final_gloscope.rds'],
                                    'rows': [{'dataset': 'Joanito',
                                              'method': 'gloscope',
                                              'parameter': '',
                                              'view': 'batch_effect_uncorrected'},
                                             {'dataset': 'Stephenson',
                                              'method': 'gloscope',
                                              'parameter': '',
                                              'view': 'batch_effect_uncorrected'},
                                             {'dataset': 'Alzheimer',
                                              'method': 'gloscope',
                                              'parameter': '',
                                              'view': 'batch_effect_uncorrected'},
                                             {'dataset': 'Breast_cancer',
                                              'method': 'gloscope',
                                              'parameter': '',
                                              'view': 'batch_effect_uncorrected'},
                                             {'dataset': 'Covid19_PBMC',
                                              'method': 'gloscope',
                                              'parameter': '',
                                              'view': 'batch_effect_uncorrected'},
                                             {'dataset': 'Kidney_KPMP_full',
                                              'method': 'gloscope',
                                              'parameter': '',
                                              'view': 'batch_effect_uncorrected'},
                                             {'dataset': 'Lupus_PBMC',
                                              'method': 'gloscope',
                                              'parameter': '',
                                              'view': 'batch_effect_uncorrected'},
                                             {'dataset': 'Lung',
                                              'method': 'gloscope',
                                              'parameter': '',
                                              'view': 'batch_effect_uncorrected'},
                                             {'dataset': 'Joanito',
                                              'method': 'mrvi',
                                              'parameter': '',
                                              'view': 'batch_effect_uncorrected'},
                                             {'dataset': 'Stephenson',
                                              'method': 'mrvi',
                                              'parameter': '',
                                              'view': 'batch_effect_uncorrected'},
                                             {'dataset': 'Alzheimer',
                                              'method': 'mrvi',
                                              'parameter': '',
                                              'view': 'batch_effect_uncorrected'},
                                             {'dataset': 'Breast_cancer',
                                              'method': 'mrvi',
                                              'parameter': '',
                                              'view': 'batch_effect_uncorrected'}],
                                    'runtime': ['embeddings/Alzheimer_batch_effect_uncorrected_final_hvg2000_highres_mrvi_dists.feather.runtime.json',
                                                'embeddings/Breast_cancer_batch_effect_uncorrected_final_hvg2000_highres_mrvi_dists.feather.runtime.json',
                                                'embeddings/Joanito_batch_effect_uncorrected_final_hvg2000_highres_mrvi_dists.feather.runtime.json',
                                                'embeddings/Stephenson_batch_effect_uncorrected_final_hvg2000_highres_mrvi_dists.feather.runtime.json'],
                                    'timings': ['execution_times_batch_effect_uncorrected_final_gloscope_Alzheimer.feather',
                                                'execution_times_batch_effect_uncorrected_final_gloscope_Breast_cancer.feather',
                                                'execution_times_batch_effect_uncorrected_final_gloscope_Covid19_PBMC.feather',
                                                'execution_times_batch_effect_uncorrected_final_gloscope_Joanito.feather',
                                                'execution_times_batch_effect_uncorrected_final_gloscope_Kidney_KPMP_full.feather',
                                                'execution_times_batch_effect_uncorrected_final_gloscope_Lung.feather',
                                                'execution_times_batch_effect_uncorrected_final_gloscope_Lupus_PBMC.feather',
                                                'execution_times_batch_effect_uncorrected_final_gloscope_Stephenson.feather',
                                                'execution_times_batch_effect_uncorrected_final_mrvi_Alzheimer.feather',
                                                'execution_times_batch_effect_uncorrected_final_mrvi_Breast_cancer.feather',
                                                'execution_times_batch_effect_uncorrected_final_mrvi_Joanito.feather',
                                                'execution_times_batch_effect_uncorrected_final_mrvi_Stephenson.feather']},
 'benchmark': {'core': ['embeddings/Bassez_hvg1000_highres_scpoli_dims15_embs.feather',
                        'embeddings/Bassez_hvg2000_highres_scpoli_dims10_embs.feather',
                        'embeddings/Bassez_hvg2000_highres_scpoli_dims2_embs.feather',
                        'embeddings/Bassez_hvg2000_highres_scpoli_dims3_embs.feather',
                        'embeddings/Bassez_hvg2000_highres_scpoli_dims5_embs.feather',
                        'embeddings/Bassez_hvg3000_highres_scpoli_dims15_embs.feather',
                        'embeddings/Stephenson_hvg1000_highres_scpoli_dims15_embs.feather',
                        'embeddings/Stephenson_hvg2000_highres_scpoli_dims10_embs.feather',
                        'embeddings/Stephenson_hvg2000_highres_scpoli_dims2_embs.feather',
                        'embeddings/Stephenson_hvg2000_highres_scpoli_dims3_embs.feather',
                        'embeddings/Stephenson_hvg2000_highres_scpoli_dims5_embs.feather',
                        'embeddings/Stephenson_hvg2000_lowres_scpoli_dims15_embs.feather',
                        'embeddings/Stephenson_hvg2000_mrvi_dists.feather',
                        'embeddings/Stephenson_hvg3000_highres_scpoli_dims15_embs.feather',
                        'results/Adams_trans.rds',
                        'results/Bassez_Pseudobulk_CT_HR_hvg2000.rds',
                        'results/Bassez_Pseudobulk_CT_HR_hvg500.rds',
                        'results/Bassez_trans.rds',
                        'results/Gongsharma_cmv_young_males_trans.rds',
                        'results/Kfoury_trans.rds',
                        'results/Kim_trans.rds',
                        'results/Lee_trans.rds',
                        'results/Pelka_trans.rds',
                        'results/Smillie_trans.rds',
                        'results/Stephenson_GloScope_hvg1000_pcadims30.rds',
                        'results/Stephenson_GloScope_hvg2000_pcadims10.rds',
                        'results/Stephenson_GloScope_hvg2000_pcadims50.rds',
                        'results/Stephenson_GloScope_hvg3000_pcadims30.rds',
                        'results/Stephenson_composition.rds',
                        'results/Stephenson_metadata.rds',
                        'results/Stephenson_mofa.rds',
                        'results/Stephenson_pseudobulk.rds',
                        'results/Stephenson_scitd.rds',
                        'results/Stephenson_trans.rds',
                        'results/Wu_trans.rds',
                        'results/Zhang_trans.rds'],
               'rows': [{'dataset': 'Stephenson', 'method': 'composition', 'parameter': '', 'view': 'benchmark_analysis'},
                        {'dataset': 'Stephenson',
                         'method': 'gloscope',
                         'parameter': 'hvg2000_pcadims10',
                         'view': 'benchmark_analysis'},
                        {'dataset': 'Stephenson',
                         'method': 'gloscope',
                         'parameter': 'hvg2000_pcadims50',
                         'view': 'benchmark_analysis'},
                        {'dataset': 'Stephenson',
                         'method': 'gloscope',
                         'parameter': 'hvg1000_pcadims30',
                         'view': 'benchmark_analysis'},
                        {'dataset': 'Stephenson',
                         'method': 'gloscope',
                         'parameter': 'hvg3000_pcadims30',
                         'view': 'benchmark_analysis'},
                        {'dataset': 'Stephenson', 'method': 'mofa', 'parameter': '', 'view': 'benchmark_analysis'},
                        {'dataset': 'Stephenson', 'method': 'mrvi', 'parameter': 'hvg2000', 'view': 'benchmark_analysis'},
                        {'dataset': 'Bassez',
                         'method': 'pseudobulk',
                         'parameter': 'CT_HR_hvg2000,CT_HR_hvg500',
                         'view': 'benchmark_analysis'},
                        {'dataset': 'Stephenson', 'method': 'pseudobulk', 'parameter': '', 'view': 'benchmark_analysis'},
                        {'dataset': 'Stephenson', 'method': 'scitd', 'parameter': '', 'view': 'benchmark_analysis'},
                        {'dataset': 'Stephenson',
                         'method': 'scpoli',
                         'parameter': 'hvg2000_lowres_dims15',
                         'view': 'benchmark_analysis'},
                        {'dataset': 'Stephenson',
                         'method': 'scpoli',
                         'parameter': 'hvg1000_highres_dims15',
                         'view': 'benchmark_analysis'},
                        {'dataset': 'Stephenson',
                         'method': 'scpoli',
                         'parameter': 'hvg2000_highres_dims2',
                         'view': 'benchmark_analysis'},
                        {'dataset': 'Stephenson',
                         'method': 'scpoli',
                         'parameter': 'hvg2000_highres_dims3',
                         'view': 'benchmark_analysis'},
                        {'dataset': 'Stephenson',
                         'method': 'scpoli',
                         'parameter': 'hvg2000_highres_dims5',
                         'view': 'benchmark_analysis'},
                        {'dataset': 'Stephenson',
                         'method': 'scpoli',
                         'parameter': 'hvg2000_highres_dims10',
                         'view': 'benchmark_analysis'},
                        {'dataset': 'Stephenson',
                         'method': 'scpoli',
                         'parameter': 'hvg3000_highres_dims15',
                         'view': 'benchmark_analysis'},
                        {'dataset': 'Bassez',
                         'method': 'scpoli',
                         'parameter': 'hvg2000_highres_dims2',
                         'view': 'benchmark_analysis'},
                        {'dataset': 'Bassez',
                         'method': 'scpoli',
                         'parameter': 'hvg2000_highres_dims3',
                         'view': 'benchmark_analysis'},
                        {'dataset': 'Bassez',
                         'method': 'scpoli',
                         'parameter': 'hvg2000_highres_dims5',
                         'view': 'benchmark_analysis'},
                        {'dataset': 'Bassez',
                         'method': 'scpoli',
                         'parameter': 'hvg2000_highres_dims10',
                         'view': 'benchmark_analysis'},
                        {'dataset': 'Bassez',
                         'method': 'scpoli',
                         'parameter': 'hvg1000_highres_dims15',
                         'view': 'benchmark_analysis'},
                        {'dataset': 'Bassez',
                         'method': 'scpoli',
                         'parameter': 'hvg3000_highres_dims15',
                         'view': 'benchmark_analysis'},
                        {'dataset': 'Stephenson', 'method': 'trans', 'parameter': '', 'view': 'benchmark_analysis'},
                        {'dataset': 'Bassez', 'method': 'trans', 'parameter': '', 'view': 'benchmark_analysis'},
                        {'dataset': 'Adams', 'method': 'trans', 'parameter': '', 'view': 'benchmark_analysis'},
                        {'dataset': 'Gongsharma_cmv_young_males',
                         'method': 'trans',
                         'parameter': '',
                         'view': 'benchmark_analysis'},
                        {'dataset': 'Kfoury', 'method': 'trans', 'parameter': '', 'view': 'benchmark_analysis'},
                        {'dataset': 'Kim', 'method': 'trans', 'parameter': '', 'view': 'benchmark_analysis'},
                        {'dataset': 'Lee', 'method': 'trans', 'parameter': '', 'view': 'benchmark_analysis'},
                        {'dataset': 'Pelka', 'method': 'trans', 'parameter': '', 'view': 'benchmark_analysis'},
                        {'dataset': 'Smillie', 'method': 'trans', 'parameter': '', 'view': 'benchmark_analysis'},
                        {'dataset': 'Wu', 'method': 'trans', 'parameter': '', 'view': 'benchmark_analysis'},
                        {'dataset': 'Zhang', 'method': 'trans', 'parameter': '', 'view': 'benchmark_analysis'}],
               'runtime': ['embeddings/Bassez_hvg1000_highres_scpoli_dims15_embs.feather.runtime.json',
                           'embeddings/Bassez_hvg2000_highres_scpoli_dims10_embs.feather.runtime.json',
                           'embeddings/Bassez_hvg2000_highres_scpoli_dims2_embs.feather.runtime.json',
                           'embeddings/Bassez_hvg2000_highres_scpoli_dims3_embs.feather.runtime.json',
                           'embeddings/Bassez_hvg2000_highres_scpoli_dims5_embs.feather.runtime.json',
                           'embeddings/Bassez_hvg3000_highres_scpoli_dims15_embs.feather.runtime.json',
                           'embeddings/Stephenson_hvg1000_highres_scpoli_dims15_embs.feather.runtime.json',
                           'embeddings/Stephenson_hvg2000_highres_scpoli_dims10_embs.feather.runtime.json',
                           'embeddings/Stephenson_hvg2000_highres_scpoli_dims2_embs.feather.runtime.json',
                           'embeddings/Stephenson_hvg2000_highres_scpoli_dims3_embs.feather.runtime.json',
                           'embeddings/Stephenson_hvg2000_highres_scpoli_dims5_embs.feather.runtime.json',
                           'embeddings/Stephenson_hvg2000_lowres_scpoli_dims15_embs.feather.runtime.json',
                           'embeddings/Stephenson_hvg2000_mrvi_dists.feather.runtime.json',
                           'embeddings/Stephenson_hvg3000_highres_scpoli_dims15_embs.feather.runtime.json',
                           'results/Adams_trans.rds.runtime.json',
                           'results/Bassez_trans.rds.runtime.json',
                           'results/Gongsharma_cmv_young_males_trans.rds.runtime.json',
                           'results/Kfoury_trans.rds.runtime.json',
                           'results/Kim_trans.rds.runtime.json',
                           'results/Lee_trans.rds.runtime.json',
                           'results/Pelka_trans.rds.runtime.json',
                           'results/Smillie_trans.rds.runtime.json',
                           'results/Stephenson_trans.rds.runtime.json',
                           'results/Wu_trans.rds.runtime.json',
                           'results/Zhang_trans.rds.runtime.json'],
               'timings': ['execution_times_composition_Stephenson.feather',
                           'execution_times_gloscope_Stephenson_hvg1000_pcadims30.feather',
                           'execution_times_gloscope_Stephenson_hvg2000_pcadims10.feather',
                           'execution_times_gloscope_Stephenson_hvg2000_pcadims50.feather',
                           'execution_times_gloscope_Stephenson_hvg3000_pcadims30.feather',
                           'execution_times_mofa_Stephenson.feather',
                           'execution_times_mrvi_Stephenson_hvg2000.feather',
                           'execution_times_pseudobulk_Bassez_selected.feather',
                           'execution_times_pseudobulk_Stephenson.feather',
                           'execution_times_scitd_Stephenson.feather',
                           'execution_times_scpoli_Bassez_hvg1000_highres_dims15.feather',
                           'execution_times_scpoli_Bassez_hvg2000_highres_dims10.feather',
                           'execution_times_scpoli_Bassez_hvg2000_highres_dims2.feather',
                           'execution_times_scpoli_Bassez_hvg2000_highres_dims3.feather',
                           'execution_times_scpoli_Bassez_hvg2000_highres_dims5.feather',
                           'execution_times_scpoli_Bassez_hvg3000_highres_dims15.feather',
                           'execution_times_scpoli_Stephenson_hvg1000_highres_dims15.feather',
                           'execution_times_scpoli_Stephenson_hvg2000_highres_dims10.feather',
                           'execution_times_scpoli_Stephenson_hvg2000_highres_dims2.feather',
                           'execution_times_scpoli_Stephenson_hvg2000_highres_dims3.feather',
                           'execution_times_scpoli_Stephenson_hvg2000_highres_dims5.feather',
                           'execution_times_scpoli_Stephenson_hvg2000_lowres_dims15.feather',
                           'execution_times_scpoli_Stephenson_hvg3000_highres_dims15.feather',
                           'execution_times_trans_Adams.feather',
                           'execution_times_trans_Bassez.feather',
                           'execution_times_trans_Gongsharma_cmv_young_males.feather',
                           'execution_times_trans_Kfoury.feather',
                           'execution_times_trans_Kim.feather',
                           'execution_times_trans_Lee.feather',
                           'execution_times_trans_Pelka.feather',
                           'execution_times_trans_Smillie.feather',
                           'execution_times_trans_Stephenson.feather',
                           'execution_times_trans_Wu.feather',
                           'execution_times_trans_Zhang.feather']}}

def md5(path):return hashlib.md5(path.read_bytes()).hexdigest()
def put(path,value):
 path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(value));return dict(path=str(path),sha256=pub.digest(path))
def fixture(base,context):
 original=base/'_ecoda_runs'/pub.MAIN_RUNS[context];run=base/'_ecoda_runs/new_release';run.mkdir(parents=True)
 paths=FIXTURE[context];root=original/'artifacts'/context
 for rel in paths['core']+paths['runtime']:
  p=root/rel;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(b'tiny scope fixture')
 for name in paths['timings']:
  p=original/'logs'/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(b'tiny timing fixture')
 catalog=dict(context=context,rows=paths['rows']);validation=dict(schema=pub.AUDIT54_SCHEMA,selected_core_paths=[str(root/r) for r in paths['core']],runtime_companions=[str(root/r) for r in paths['runtime']],timing_paths=[str(original/'logs'/r) for r in paths['timings']])
 return original,run,catalog,validation

class Tests(unittest.TestCase):
 def test_exact54_public_scope_and147_companions(self):
  with tempfile.TemporaryDirectory() as tmp:
   base=pathlib.Path(tmp).resolve();totals=[0,0,0,0]
   for context in FIXTURE:
    original,run,catalog,validation=fixture(base/context.replace('/','_'),context);core,runtime,timing=pub.audit54_paths(validation,catalog,context,original)
    observed=(len(catalog['rows']),len(core),len(runtime),len(timing));self.assertEqual(observed,pub.AUDIT54_COUNTS[context]);totals=[a+b for a,b in zip(totals,observed)]
    bad=copy.deepcopy(catalog);bad['rows'][0]['method']='zeroimp'
    with self.assertRaises(ValueError):pub.audit54_paths(validation,bad,context,original)
    bad=copy.deepcopy(validation);bad['runtime_companions'].pop()
    with self.assertRaises(ValueError):pub.audit54_paths(bad,catalog,context,original)
   self.assertEqual(totals,[54,56,37,54])
 def test_filtered_scitd14_and_null_own_labels(self):
  ids=['s'+str(i) for i in range(53)];truth={x:'A' if i%2 else 'B' for i,x in enumerate(ids)};drop=ids[:14];retained=ids[14:]
  def bundle(samples,dropped):return dict(sample_ids=samples,dropped_sample_ids=dropped,stored_labels=[truth[x] for x in samples],timing=dict(method='scITD_hvg2000_factors5',time_secs=1.))
  observed=dict(kind='result_bundles',bundles={'scITD_hvg2000_factors5':bundle(retained,drop)})
  self.assertEqual(pub.audit54_bundle_ids(observed,ids,truth),(retained,drop))
  bad=copy.deepcopy(observed);bad['bundles']['scITD_hvg2000_factors5']['sample_ids'].reverse()
  with self.assertRaises(ValueError):pub.audit54_bundle_ids(bad,ids,truth)
  b=bundle(ids,[]);b['stored_labels'].reverse();b['timing']['method']='ECODA_authors_HR_NULL'
  self.assertEqual(pub.audit54_bundle_ids(dict(kind='result_bundles',bundles={'ECODA_authors_HR_NULL':b}),ids,truth),(ids,[]))
  b['stored_labels'][0]='FOREIGN'
  with self.assertRaises(ValueError):pub.audit54_bundle_ids(dict(kind='result_bundles',bundles={'ECODA_authors_HR_NULL':b}),ids,truth)
 def test_transformation_named_package_is_not_faked_positional_order(self):
  ids=['s1','s2','s3'];truth={'s1':'A','s2':'B','s3':'A'}
  observed=dict(kind='trans_v2',policy={'version':2},draws=120,n_samples=3,stored_labels=['B','A','A'],disposition='EXPLICIT_TRANS_V2_POLICY_VALIDATED_PENDING_ROOT_REVIEW')
  self.assertEqual(pub.audit54_bundle_ids(observed,ids,truth),(ids,[]))
  observed['disposition']='SUMMARY_ONLY_NO_STORED_LABELS_OR_DISTANCES_NEEDS_SEPARATE_POLICY_REVIEW'
  with self.assertRaises(ValueError):pub.audit54_bundle_ids(observed,ids,truth)
 def test_observed_stat_proof_and_absent_owner_are_immutable(self):
  with tempfile.TemporaryDirectory() as tmp:
   base=pathlib.Path(tmp).resolve();p=base/'tiny.runtime.json';p.write_text('{}');side=base/'proof';side.write_text('proof');record=base/'record';record.write_text('record');owner=base/'absent_owner'
   stat=lambda p:{k:getattr(p.stat(),k) for k in (*pub.STAT_KEYS,'st_mode')}
   row=dict(path=str(p),bytes_valid=True,proofs_valid=True,stat_after=stat(p),sha256=pub.digest(p),sidecar=dict(path=str(side),sha256=pub.digest(side),stat=stat(side)),producer_record=dict(path=str(record),sha256=pub.digest(record),stat=stat(record)),owner=dict(path=str(owner),exists=False))
   pub.audit54_observed_file(row);owner.write_text('appeared')
   with self.assertRaises(ValueError):pub.audit54_observed_file(row)
   owner.unlink();p.write_text('changed')
   with self.assertRaises(ValueError):pub.audit54_observed_file(row)
 def test_copied_corrected_core_original_runtime_usable_by_current_P6(self):
  import pandas as pd
  with tempfile.TemporaryDirectory() as tmp:
   base=pathlib.Path(tmp).resolve();context='batch_effect/corrected_final';original,run,catalog,validation=fixture(base,context);root=original/'artifacts'/context
   runtime=pathlib.Path(validation['runtime_companions'][0]);core=pathlib.Path(str(runtime)[:-13]);ids=['s1','s2'];pd.DataFrame({'Sample':ids,'s1':[0.,1.],'s2':[1.,0.]}).to_feather(core)
   value=dict(schema_version=1,artifact_path=str(core),artifact_md5=md5(core),dataset='Alzheimer',method='MrVI_hvg2000',time_secs=2.,mem_GB=None,batch_contract={'native_original_identity':'fixture'})
   runtime.write_text(json.dumps(value));before={str(p):p.read_bytes() for p in (core,runtime)}
   cref=put(base/'catalog.json',catalog);vref=put(base/'validation.json',validation);science=put(base/'science.json',dict(catalog=cref,scientific_validation=vref))
   control=run/'manifests/control.json';put(control,dict(schema=pub.REBOOT_SCHEMA,run_id=run.name,context=context,science_receipt=science));dest=run/'artifacts'/context;nas=base/'nas';nas.mkdir();mappings=[];inventory=[]
   for source in (core,runtime):
    relative=source.relative_to(root);targets=[]
    for boundary,target in ((run,dest/relative),(nas,nas/relative)):
     target.parent.mkdir(parents=True,exist_ok=True);pub.copy_payload(SimpleNamespace(source=str(source),destination=str(target),boundary=str(boundary),md5=md5(source),size=str(source.stat().st_size),control=str(control),control_md5=md5(control)));targets.append(str(target))
    kind='runtime' if source==runtime else 'result';item=dict(source=str(source),kind=kind,md5=md5(source));inventory.append(item);mappings.append(dict(item,scratch=targets[0],nas=targets[1]))
   lineage=pub.audit54_runtime_lineage(inventory,mappings);self.assertEqual(len(lineage),1);self.assertEqual(json.loads(pathlib.Path(mappings[1]['scratch']).read_text())['artifact_path'],str(core));self.assertEqual(before,{str(p):p.read_bytes() for p in (core,runtime)})
   # Producer-path validation must NOT be claimed for byte-identical copies.
   native_spec=importlib.util.spec_from_file_location('native_contract',ROOT/'src/5_run_benchmark_methods/matrix_artifact_validator.py');native=importlib.util.module_from_spec(native_spec);native_spec.loader.exec_module(native)
   copied=pathlib.Path(mappings[0]['scratch']);copied_runtime=pathlib.Path(mappings[1]['scratch'])
   for p in (copied,copied_runtime):pathlib.Path(str(p)+'.md5').write_text('MD5='+md5(p)+'\nSIZE='+str(p.stat().st_size)+'\nPATH='+str(p)+'\n')
   with self.assertRaisesRegex(ValueError,'wrong artifact_path'):native._read_feather_batch_contract(copied,require_summary=False)
   evidence=base/'checked_lineage.json';evidence.write_text(json.dumps(lineage));script=base/'p6.R';script.write_text(R_P6)
   env=dict(os.environ,R_HOME=str(ROOT/'.pixi/envs/default/lib/R'));result=subprocess.run([str(ROOT/'.pixi/envs/default/bin/Rscript'),'--vanilla',str(script),str(ROOT),str(base),str(copied),str(evidence)],env=env,capture_output=True,text=True)
   self.assertEqual(result.returncode,0,result.stdout+result.stderr);self.assertIn('CURRENT_P6_READER_PROJECTION_OK',result.stdout)
   bad=base/'arbitrary.runtime.json';bad.write_text('{}')
   with self.assertRaises(ValueError):pub.copy_payload(SimpleNamespace(source=str(bad),destination=str(nas/'bad.runtime.json'),boundary=str(nas),md5=md5(bad),size=str(bad.stat().st_size),control=str(control),control_md5=md5(control)))
R_P6=r"""
a<-commandArgs(TRUE);source(file.path(a[1],'src/utils/analysis_release.R'));root<-normalizePath(a[2]);md5<-function(p)unname(tools::md5sum(p))
row<-data.frame(dataset='Alzheimer',view='batch_effect_corrected',method_key='MrVI_hvg2000',kind='distance_feather',artifact_path=a[3],bundle_key='',run_root=dirname(dirname(dirname(a[3]))),md5=md5(a[3]),receipt_path=file.path(root,'receipt.rds'),receipt_md5='',sample_id_col='Sample',stringsAsFactors=FALSE)
# This synthetic review is of the checked publication lineage, not a native new-fit claim.
receipt<-c(as.list(row[c('dataset','view','method_key','kind','artifact_path','bundle_key','run_root','md5')]),list(schema='ecoda_analysis_release_receipt_v1',status='reviewed',origin='legacy_reuse',input_md5=strrep('a',32),source_policy='reviewed_original_fit_copy',evidence='synthetic checked original/runtime/copy lineage',evidence_path=a[4],evidence_md5=md5(a[4])))
saveRDS(receipt,row$receipt_path);row$receipt_md5<-md5(row$receipt_path)
scope<-row[c('dataset','view','method_key')];scope$require_new<-FALSE;scope$current_input_md5<-strrep('a',32);scope$source_policy<-'reviewed_original_fit_copy';scope$status<-'approved'
control<-list(ledger=row,expected_scope=scope);control_path<-file.path(root,'control.rds');saveRDS(control,control_path);s<-file.path(root,'source');r<-file.path(root,'runtime');writeLines('fixture source',s);writeLines('fixture runtime',r)
value<-load_analysis_release(row,scope,require_projection=TRUE)[[1]]$result;stopifnot(identical(rownames(as.matrix(value$dist_mat)),c('s1','s2')))
projected<-project_analysis_release(control,file.path(root,'thin'),c(source=s,runtime=r,control=control_path));loaded<-load_analysis_release(projected$ledger,projected$expected_scope,require_projection=TRUE)[[1]]
stopifnot(identical(loaded$receipt$origin,'legacy_reuse'),identical(as.matrix(loaded$result$dist_mat),as.matrix(value$dist_mat)))
bad<-projected$expected_scope;bad$require_new<-TRUE;stopifnot(inherits(tryCatch(load_analysis_release(projected$ledger,bad,require_projection=TRUE),error=identity),'error'))
cat('CURRENT_P6_READER_PROJECTION_OK\n')
"""
if __name__=='__main__':unittest.main()
