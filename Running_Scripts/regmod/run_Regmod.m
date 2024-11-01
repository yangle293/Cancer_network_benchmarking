function run_Regmod(n)

n = n+1;
data_path = './data/';
network_path = './network/';
result_path = './result/';
data_version_list = {'corum_0','corum_1','corum_2','corum_3','corum_4',...
    'reactome_0','reactome_1','reactome_2','reactome_3','reactome_4'};
network_list = {'biogrid','irefindex18','reactome21','string'};
networks = repelem(network_list,1,1100);

score = allcomb(data_version_list,num2cell(1:10),num2cell(0.01:0.01:0.11));
scores = repmat(score,length(network_list),1);

% indexing
s = [scores{n,1},'_z_beta_',num2str(scores{n,3}),'_',num2str(scores{n,2})];
network = networks{n};

s
network
%% load data
[index_start,index_end] = import_edge_list([network_path,network,'_edge_list']);
[index,gene] = import_gene_index([network_path,network,'_index_gene']);
[Gene_score,Z] = import_gene_score([data_path,s,'.txt']);

%% precalculate
G = graph(gene(index_start),gene(index_end));
d = full(G.adjacency);
deg = sum(d); 
[C,IA,IB] = intersect(G.Nodes.Name,Gene_score,'stable');
z = min(Z)*ones(length(G.Nodes.Name),1);
z(IA) = Z(IB);


beta = [1,2,3,4,5];
C_list = [1,5,10];
f_list = [2,4,6];
dfm = cell(length(beta),1);
Lapd = -1*diag(deg) + d;
for i = 1:length(beta)
    dfm{i} = padm(beta(i)*Lapd);
end


%result_f = cell(length(beta),length(C_list));
for i = 1:length(beta)
    for j = 1:length(C_list)
        z_smooth = RegMOD(d, z, dfm{i}, C_list(j));
        m = mean(z_smooth);ss = mad(z_smooth);
        for k = 1:length(f_list)
            f = f_list(k);
            ind = find(z_smooth>(m+f*ss));
            if isempty(ind)
                continue
            end
            tmp_h = G.subgraph(ind);cc = conncomp(tmp_h);t = tabulate(cc);
            v = t(find(t(:,2)>1),1);
            if isempty(v)
                continue
            end
            indx = [];
            networks = cell(length(v),1);
            for l = 1:length(v)
                networks{l} = tmp_h.Nodes.Name(find(cc==v(l)));
            end
            net = cell2table(cellfun(@strjoin,networks,'UniformOutput',0));
            filename = [result_path,network,s,'_b',num2str(beta(i)),'_c',num2str(C_list(j)),'_f',num2str(f_list(k)),'.txt']
            writetable(net,filename,'WriteVariableNames',false);
        end
    end
end

%save([result_path,network,s,'_result.mat'],"result_f")

end

