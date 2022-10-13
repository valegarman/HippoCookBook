function mtint = checkPosScaling(mtint)
% checks the pos file header and the max and min x/y values to make sure
% they will fit in. if not then the positional data is scaled appropriately

mtint.header;