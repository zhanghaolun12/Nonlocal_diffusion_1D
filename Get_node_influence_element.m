function [element_index] = Get_node_influence_element(node_index,T_partition)
%ȷ������Ԫ�ڵ��Ӱ�쵥Ԫ������������֧��
element_index = []; %Ӱ�쵥Ԫָ��
number_of_elements = size(T_partition, 2); %��Ԫ��
for n = 1 : number_of_elements
    if ~isempty(find(T_partition(:,n)==node_index))
        element_index = [element_index, n];
    end
end

end