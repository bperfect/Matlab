function result = convert2rho(input,format)
%% convert a field variable in psi coordinates to rho coordinates
%% should work as long as the 2-leading dimensions are x and y, respectively

nd = ndims(input);
if nd > 2
    index = repmat({':'},1,nd-2);
else
    index = {};
end

dims = size(input);

%%case psi
if format == 'psi'
result = input(2:end,2:end,index{:});
for row = 1:dims(1)-1
    for col = 1:dims(2)-1
        result(row,col,index{:}) = 1/4*(input(row,col,index{:})...
            +input(row+1,col,index{:})+input(row,col+1,index{:})...
            +input(row+1,col+1,index{:}));
    end
end
end

%%case u
if format=='u'
result = input(2:end,2:end-1,index{:});
for row = 1:dims(1)-1
    result(row,:,index{:}) = 1/2*(input(row+1,2:end-1,index{:})+input(row,2:end-1,index{:}));
end
end

%%case u
if format=='v'
result = input(2:end-1,2:end,index{:});
for col = 1:dims(2)-1
    result(:,col,index{:}) = 1/2*(input(2:end-1,col+1,index{:})+input(2:end-1,col,index{:}));
end
end

%case rho
if format == 'rho'
    result = input(2:end-1,2:end-1,index{:});
end



