tic
bigMatrixA = rand(3000000,80);
bigMatrixB = rand(80,30);
bigMatrixC = bigMatrixA * bigMatrixB;
toc
disp("done");
