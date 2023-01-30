#!/usr/bin/env python 
import Data_DB_manager as db_mng
import Eventalign_DB as db
import TranscriptObject
import Sample
import TxComp

sample2condition = {'KO1':"KO", 'KO2':"KO", 'KO3':"KO", 
                    'WT1':"WT", 'WT2':"WT", 'WT3':"WT"}
#test dataset
test_samples = {'KO1':'/hps/nobackup/birney/users/logan/nanocompore_rewrite/sqlite_KO1_eventalign_collapse.db',
                'WT1':'/hps/nobackup/birney/users/logan/nanocompore_rewrite/WT1_eventalign_collapse.db'}

test_reference_samples = {'KO1':'/hps/nobackup/birney/users/logan/nanocompore_rewrite/sqlite_KO1_eventalign_collapse.db',
                          'KO2':'/hps/nobackup/birney/users/logan/nanocompore_rewrite/sqlite_KO1_eventalign_collapse2.db'}

test_test_samples = {'WT1':'/hps/nobackup/birney/users/logan/nanocompore_rewrite/WT1_eventalign_collapse.db',
                     'WT2':'/hps/nobackup/birney/users/logan/nanocompore_rewrite/WT1_eventalign_collapse2.db'}

#list of transcripts for testing purposes
tx = ['YCR010C_mRNA', 'YNL208W_mRNA', 'YNL037C_mRNA',
      'YDR099W_mRNA', 'YMR303C_mRNA', 'YNL155W_mRNA',
      'YGL009C_mRNA', 'YEL060C_mRNA', 'YBR067C_mRNA',
      'YPR036W-A_mRNA', 'YCR005C_mRNA']
test_tx_name = tx[0]
tx_length = 852
max_coverage = float('inf')
min_coverage = 30
pos_0 = 1

def test_db_manager():
      #test that db is properly opened
    for label in test_samples:
        test_db_path = test_samples[label]
        test_db = db_mng.Data_DB_manager(test_db_path)

        #test that data for a particular trasncript can be gathered
        test_data = test_db.getData(test_tx_name)

        #test that the correct number of positions and reads for a particular position are correct
        print(label)
        print('number of positions', len(test_data.keys()))
        print('dwell_times', len(test_data[pos_0]['dwell_times']))
        print('intensities', len(test_data[pos_0]['intensities']))

        #test that the db can be closed
        test_db.closeDB()
        print('\n')


def test_sample_object():
    #test that a sample object is created properly
    for sample_label in test_samples:
        cond_label = sample2condition[sample_label]
        db_path = test_samples[sample_label]

        sample = Sample.Sample(cond_label, sample_label, db_path, tx[0], max_coverage)
        valid_positions = sample.getValidPositions(min_coverage)
        print('valid positions', len(valid_positions))
        if pos_0 in valid_positions:
            print(f"pos {pos_0}")
            print('read_ids', len(sample.getReadIds(pos_0)), sample.getReadIds(pos_0)[0:10])
            print('intensities', len(sample.getKmerIntensityData(pos_0)), sample.getKmerIntensityData(pos_0)[0:10])
            print('dwell_times', len(sample.getKmerDwelldata(pos_0)), sample.getKmerDwelldata(pos_0)[0:10])
            print('sample_ids', len(sample.getKmerLabelData(pos_0)), sample.getKmerLabelData(pos_0)[0:10])

        #test that a sample object closes the db properly
        sample.closeDB()
        print('\n')


def test_transcript_object():
    #test that a transcript object functions as expected
    test_transcript_data = TranscriptObject.Transcript_Data(test_tx_name, test_reference_samples, test_test_samples, tx_length, max_coverage, min_coverage)

    print("transcript object public variables")
    print(f"transcirpt name: {test_transcript_data.name}")
    print(f"transcript length: {test_transcript_data.length}")
    print("transcript object private variables")
    print(f"reference samples: \n{test_transcript_data._reference_samples}")
    print(f"test samples: \n{test_transcript_data._test_samples}")
    print(f"all sample labels: \n{test_transcript_data._all_sample_labels}")
    print(f"min coverage: {test_transcript_data._min_coverage}")
    for coverage_level in [min_coverage, 50, 10000]:
        if test_transcript_data.enoughTxCoverage(coverage_level):
            print(f'Greater than {coverage_level} reads')
        else:
            print(f'Less than {coverage_level} reads')

        print("number of valid positions", len(test_transcript_data.getValidPositions()))
        print("Reference intensity data")
        print(len(test_transcript_data.getReferenceIntensityData(pos_0)))
        print(type(test_transcript_data.getReferenceIntensityData(pos_0)))
        print("Test intensity data")
        print(len(test_transcript_data.getTestIntensityData(pos_0)))
        print(type(test_transcript_data.getTestIntensityData(pos_0)))

        print("Reference dwell data")
        print(len(test_transcript_data.getReferenceDwellData(pos_0)))
        print(type(test_transcript_data.getReferenceDwellData(pos_0)))
        print("Test dwell data")
        print(len(test_transcript_data.getTestDwellData(pos_0)))
        print(type(test_transcript_data.getTestDwellData(pos_0)))

        print("Reference label data")
        print(len(test_transcript_data.getReferenceLabelData(pos_0)))
        print(type(test_transcript_data.getReferenceLabelData(pos_0)))
        print("Test label data")
        print(len(test_transcript_data.getTestLabelData(pos_0)))
        print(type(test_transcript_data.getTestLabelData(pos_0)))

        print("Reference Read Ids")
        print(len(test_transcript_data.getReferenceReadIdData(pos_0)))
        print(type(test_transcript_data.getReferenceReadIdData(pos_0)))
        print("Test Read Ids")
        print(len(test_transcript_data.getTestReadIdData(pos_0)))
        print(type(test_transcript_data.getTestReadIdData(pos_0)))


    test_transcript_data.closeAllDbs()


def test_txComp():
    #test TxComp
    txComp_object = TxComp.TxComp(num_reference_samples = 1,
                                  num_test_samples = 1, 
                                  methods=['KS', 'MW', 'TT'])
    print('Made a txcomp object')
    test_transcript_data = TranscriptObject.Transcript_Data(test_tx_name, test_reference_samples, test_test_samples, tx_length, max_coverage, min_coverage)
    txComp_object.txCompare(test_transcript_data)



    test_transcript_data.closeAllDbs()


#test_db_manager()
#test_sample_object()
#test_transcript_object()
test_txComp()