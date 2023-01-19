#!/usr/bin/env python 
import Data_DB_manager as db_mng
import Eventalign_DB as db
import TranscriptObject
import Sample

#test dataset
test_db_path = '/hps/nobackup/birney/users/logan/nanocompore_rewrite/sqlite_KO1_eventalign_collapse.db'

#test that db is properly opened
test_db = db_mng.Data_DB_manager(test_db_path)

#list of transcripts for testing purposes
tx = ['YCR010C_mRNA', 'YNL208W_mRNA', 'YNL037C_mRNA',
      'YDR099W_mRNA', 'YMR303C_mRNA', 'YNL155W_mRNA',
      'YGL009C_mRNA', 'YEL060C_mRNA', 'YBR067C_mRNA',
      'YPR036W-A_mRNA', 'YCR005C_mRNA']

#test that data for a particular trasncript can be gathered
test_data = test_db.getData(tx[0])

#test that the correct number of positions and reads for a particular position are correct
print('number of positions', len(test_data.keys()))
pos_0 = 1
print('dwell_times', len(test_data[pos_0]['dwell_times']))
print('intensities', len(test_data[pos_0]['intensities']))

#test that the db can be closed
test_db.closeDB()

#test that a sample object is created properly
test_cond_label = 'KO'
test_sample_label = 'KO1'
max_coverage = 1000
min_coverage = 30

test_sample = Sample.Sample(test_cond_label, test_sample_label, test_db_path, tx[0], max_coverage)
valid_positions = test_sample.getValidPositions(min_coverage)
print('valid positions', len(valid_positions))
if pos_0 in valid_positions:
      print('read_ids', len(test_sample.getReadIds(pos_0)), test_sample.getReadIds(pos_0)[0:10])
      print('intensities', len(test_sample.getKmerIntensityData(pos_0)), test_sample.getKmerIntensityData(pos_0)[0:10])
      print('dwell_times', len(test_sample.getKmerDwelldata(pos_0)), test_sample.getKmerDwelldata(pos_0)[0:10])
      print('sample_ids', len(test_sample.getKmerLabelData(pos_0)), test_sample.getKmerLabelData(pos_0)[0:10])

#test that a sample object closes the db properly
test_sample.closeDB()

#test that a transcript object functions as expected
test_tx_name = tx[0]
test_reference_samples = {'KO1':test_db_path}
test_test_samples = {'WT1':'/hps/nobackup/birney/users/logan/nanocompore_rewrite/WT1_eventalign_collapse.db'}
tx_length = 852
test_transcript_data = TranscriptObject.Transcript_Data(test_tx_name, test_reference_samples, test_test_samples, tx_length, max_coverage, min_coverage)

print(test_transcript_data.name, test_transcript_data.length)
if test_transcript_data.enoughTxCoverage():
      print('enough coverage')
else:
      print('not enough coverage')

if test_transcript_data.enoughTxCoverage(10000):
      print('more than 10000 reads')
else:
      print('less than 10000 reads')

print("number of valid positions", len(test_transcript_data.getValidPositions()))


test_transcript_data.closeAllDbs()