# Gaze Contingent Avoidance Study 2: Attentional Competition Version

nReps(learningtrials) = 2, for Testing: Selected rows: randchoice(16, 4)

nReps(trials) = 2, for Testing: randchoice(16, 4)

nReps(blocks) = 1

nReps(testtrials) = 2, for Testing: randchoice(24, 4)

nReps(testtrials_novelty) = 1
Selected rows: np.concatenate((np.random.permutation(np.arange(0, 24))[0:8], np.random.permutation(np.arange(24, 48))[0:10], np.random.permutation(np.arange(48, 72))[0:10], np.random.permutation(np.arange(72, 96))[0:10], np.random.permutation(np.arange(96, 120))[0:10]))
for Testing: randchoice(120, 4)