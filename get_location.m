function [ element_index ] = get_location( point, P, T )
%GET_LOCATION �ҵ�point��Ӧλ�õĵ�Ԫ�ı��
number_of_element = size(T,2);
element_index =-1;
%ͨ���������е�Ԫ��ȷ���õ����ڵ�Ԫ�ı��
for n = 1: number_of_element     
    vertices = P(:, T(:, n));   %��i����Ԫ�����Ҷ˵�����
    if(point>=vertices(1) && point<=vertices(2))
        element_index = n;
    end
end

end

