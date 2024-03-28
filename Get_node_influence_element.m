function [element_index] = Get_node_influence_element(node_index,T_partition)
%确定有限元节点的影响单元，即基函数的支集
element_index = []; %影响单元指标
number_of_elements = size(T_partition, 2); %单元数
for n = 1 : number_of_elements
    if ~isempty(find(T_partition(:,n)==node_index))
        element_index = [element_index, n];
    end
end

end