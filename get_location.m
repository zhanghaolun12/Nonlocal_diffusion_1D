function [ element_index ] = get_location( point, P, T )
%GET_LOCATION 找到point对应位置的单元的编号
number_of_element = size(T,2);
element_index =-1;
%通过遍历所有单元来确定该点所在单元的编号
for n = 1: number_of_element     
    vertices = P(:, T(:, n));   %第i个单元的左右端点坐标
    if(point>=vertices(1) && point<=vertices(2))
        element_index = n;
    end
end

end

