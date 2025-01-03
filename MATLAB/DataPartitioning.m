% Load the dataset from the CSV file
C = readmatrix('uni_167_Tomic-Stoja_RMT_no-ref-ch.ods');

% Define the sizes of each segment
sizes = [2500, 2500, 12249, 2500, 2500];

% Define the start index for each segment
start_index = [1, 2501, 5001, 17250, 19750];

% Split the dataset into segments
r1 = C(start_index(1):start_index(1)+sizes(1)-1, :);
r2 = C(start_index(2):start_index(2)+sizes(2)-1, :);
r3 = C(start_index(3):start_index(3)+sizes(3)-1, :);
r4 = C(start_index(4):start_index(4)+sizes(4)-1, :);
r5 = C(start_index(5):start_index(5)+sizes(5)-1, :);

save('Seizure_Timesig_Sub_uni20.mat', 'C', 'r1', 'r2', 'r3', 'r4', 'r5');