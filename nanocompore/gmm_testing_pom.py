#!/usr/bin/env python 
import Data_DB_manager as db_mng
import Eventalign_DB as db
import TranscriptObject
import Sample
import TxComp

import sklearn
import sklearn.datasets 

import matplotlib.pyplot as plt
import seaborn; seaborn.set_style('whitegrid')
import numpy as np

from pomegranate import *

numpy.random.seed(0)
numpy.set_printoptions(suppress=True)


sample2condition = {'KO1':"KO", 'KO2':"KO", 'KO3':"KO", 
                    'WT1':"WT", 'WT2':"WT", 'WT3':"WT"}

#test dataset
test_samples = {'KO1':'/hps/nobackup/birney/users/logan/nanocompore_rewrite/sqlite_KO1_eventalign_collapse.db',
                'WT1':'/hps/nobackup/birney/users/logan/nanocompore_rewrite/WT1_eventalign_collapse.db'}

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

def one_sample_per_condition():
    test_reference_samples = {'KO1':'/hps/nobackup/birney/users/logan/nanocompore_rewrite/sqlite_KO1_eventalign_collapse.db'}

    test_test_samples = {'WT1':'/hps/nobackup/birney/users/logan/nanocompore_rewrite/WT1_eventalign_collapse.db'}

    return test_reference_samples, test_test_samples

def two_sample_per_condition():
    test_reference_samples = {'KO1':'/hps/nobackup/birney/users/logan/nanocompore_rewrite/sqlite_KO1_eventalign_collapse.db',
                              'KO2':'/hps/nobackup/birney/users/logan/nanocompore_rewrite/sqlite_KO1_eventalign_collapse2.db'}

    test_test_samples = {'WT1':'/hps/nobackup/birney/users/logan/nanocompore_rewrite/WT1_eventalign_collapse.db',
                         'WT2':'/hps/nobackup/birney/users/logan/nanocompore_rewrite/WT1_eventalign_collapse2.db'}

    return test_reference_samples, test_test_samples


def test_transcript_object():
    #test that a transcript object functions as expected
    test_reference_samples, test_test_samples = one_sample_per_condition()
    #test_reference_samples, test_test_samples = two_sample_per_condition()

    transcript = TranscriptObject.Transcript_Data(test_tx_name, test_reference_samples, test_test_samples, tx_length, max_coverage, min_coverage)

    intensity = np.concatenate([transcript.getReferenceIntensityData(pos_0)] + [transcript.getTestIntensityData(pos_0)])
    dwell = np.concatenate([transcript.getReferenceDwellData(pos_0)] + [transcript.getTestDwellData(pos_0)])
    dwell = np.log10(dwell)

    X = np.vstack((intensity, dwell)).T 

    print(X)

    plt.figure(figsize=(8, 6))
    plt.scatter(X[:,0], X[:,1])
    plt.savefig('test_pom_plots.png')
    plt.close()

    model = GeneralMixtureModel.from_samples(MultivariateGaussianDistribution, 2, X)

    x = numpy.arange(70, 130, 0.5)
    y = numpy.arange(-3, 0, 0.025)

    xx, yy = numpy.meshgrid(x, y)
    print(xx)
    print()
    print(yy)
    print()
    x_ = numpy.array(list(zip(xx.flatten(), yy.flatten())))

    print(x_)
    p1 = MultivariateGaussianDistribution.from_samples(X).probability(x_).reshape(len(x), len(y))
    p2 = model.probability(x_).reshape(len(x), len(y))


    plt.figure(figsize=(14, 6))
    plt.subplot(121)
    plt.contourf(xx, yy, p1, cmap='Blues', alpha=0.8)
    plt.scatter(X[:,0], X[:,1])

    plt.subplot(122)
    plt.contourf(xx, yy, p2, cmap='Blues', alpha=0.8)
    plt.scatter(X[:,0], X[:,1])
    plt.savefig('test_pom_plots_2_components.png')


    transcript.closeAllDbs()



'''
#-------------------------------------------------------------------------------
#Get data from somewhere (moons data is nice for examples)
Xmoon, ymoon = sklearn.datasets.make_moons(200, shuffle = False, noise=.05, random_state=0)
Moon1 = Xmoon[:100] 
Moon2 = Xmoon[100:] 
MoonsDataSet = Xmoon

#Weight the data from moon2 much higher than moon1:
MoonWeights = numpy.array([numpy.ones(100), numpy.ones(100)*10]).flatten()

print(Moon1)
print(Moon2)
print(MoonWeights)

#Make the GMM model using pomegranate
model = pomegranate.gmm.GeneralMixtureModel.from_samples(
    pomegranate.MultivariateGaussianDistribution,   #Either single function, or list of functions
    n_components=6,     #Required if single function passed as first arg
    X=MoonsDataSet,     #data format: each row is a point-coordinate, each column is a dimension
    )

#Force the model to train again, using additional fitting parameters
model.fit(
    X=MoonsDataSet,         #data format: each row is a coordinate, each column is a dimension
    weights = MoonWeights,  #List of weights. One for each point-coordinate
    stop_threshold = .001,  #Lower this value to get better fit but take longer. 
                            #   (sklearn likes better/slower fits than pomegrante by default)
    )

'''
test_transcript_object()