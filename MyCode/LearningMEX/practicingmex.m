inputArg = int16(A);
if ~strcmp(class(inputArg),'double')
    inputArg = double(inputArg);
end
B = arrayProduct(s,inputArg)