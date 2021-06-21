function BondEnergyCom(directory, Name, Atom_info, Bond_info, Atom1, Atom2, Type, Charged_ID)
    % Create two txt files of two fragments for gaussian bond energy calculation.
    % Only limited to linear molecule (no ring)
    % Type: 1 for homolysis; 2 for heterolysis
    % Charged_ID: Protonated atom index (optional)
    
    if ~exist('Charged_ID','var')
        ChargedID = 0;
    end
    
    % Create linkage matrix to identify atoms in each fragment
    Linkage = zeros(size(Atom_info,2)); % Linkage index
    
    for i = 1:size(Bond_info,2)
       Linkage(Bond_info(i).Atom1,Bond_info(i).Atom2)=1;
       Linkage(Bond_info(i).Atom2,Bond_info(i).Atom1)=1;
    end
    
    Attached_list1 = [Atom1];
    Idle_list=[Atom1];
    
    while ~isempty(Idle_list) 
        target = Idle_list(1);
        a = find(Linkage(:,target)==1);
        for i = 1:length(a)
           if a(i)~=Atom2 && ~any(Attached_list1==a(i))
               Idle_list = [Idle_list, a(i)];
               Attached_list1 = [Attached_list1, a(i)];
           end
        end
        Idle_list(1)=[];
    end
    Attached_list1 = sort(Attached_list1);
    Attached_list2 = setdiff([1:size(Atom_info,2)], Attached_list1);
    
    % Determine charge and multiplicity
    keySet = {'H','C','N','O','F','P','S','Cl'};
    valueSet = {1,6,7,8,9,15,16,17};
    ElectronMap = containers.Map(keySet,valueSet);
    
    switch Type
        case 1
            filename1 = [Name,'_Bond',num2str(Atom1),'_',num2str(Atom2),'_Homolysis_',num2str(Atom1),'part'];
            filename2 = [Name,'_Bond',num2str(Atom1),'_',num2str(Atom2),'_Homolysis_',num2str(Atom2),'part'];
            total1=0;   total2=0;   % Number of electron
            Charge1=0;  Charge2=0;

            for i = 1:length(Attached_list1)
                total1 = total1+ElectronMap(Atom_info(Attached_list1(i)).element);
                if Attached_list1(i)==Charged_ID
                   total1 = total1-1;   % Protonation site
                   Charge1=Charge1+1;
                end
            end
            for i = 1:length(Attached_list2)
                total2 = total2+ElectronMap(Atom_info(Attached_list2(i)).element);
                if Attached_list2(i)==Charged_ID
                   total2 = total2-1;   % Protonation site
                   Charge2=Charge2+1;
                end
            end    
            Mul1 = rem(total1,2)+1;
            Mul2 = rem(total2,2)+1;
            WritePartCom(directory, [filename1,'.com'], Atom_info, Attached_list1, Charge1, Mul1);
            WritePartCom(directory, [filename2,'.com'], Atom_info, Attached_list2, Charge2, Mul2);
            
        case 2
            filename1 = [Name,'_Bond',num2str(Atom1),'_',num2str(Atom2),'_Heterolysis_to_',num2str(Atom1),'_',num2str(Atom1),'part'];
            filename2 = [Name,'_Bond',num2str(Atom1),'_',num2str(Atom2),'_Heterolysis_to_',num2str(Atom1),'_',num2str(Atom2),'part'];
            filename3 = [Name,'_Bond',num2str(Atom1),'_',num2str(Atom2),'_Heterolysis_to_',num2str(Atom2),'_',num2str(Atom1),'part'];
            filename4 = [Name,'_Bond',num2str(Atom1),'_',num2str(Atom2),'_Heterolysis_to_',num2str(Atom2),'_',num2str(Atom2),'part'];
            
            total1=1;   total2=-1;   total3=-1;   total4=1;   % Number of electron
            Charge1=-1;  Charge2=1;  Charge3=1;  Charge4=-1;

            for i = 1:length(Attached_list1)
                total1 = total1+ElectronMap(Atom_info(Attached_list1(i)).element);
                total3 = total3+ElectronMap(Atom_info(Attached_list1(i)).element);
                if Attached_list1(i)==Charged_ID
                   total1 = total1-1;   % Protonation site
                   total3 = total3-1;
                   Charge1=Charge1+1;
                   Charge3=Charge3+1;
                end
            end
            for i = 1:length(Attached_list2)
                total2 = total2+ElectronMap(Atom_info(Attached_list2(i)).element);
                total4 = total4+ElectronMap(Atom_info(Attached_list2(i)).element);
                if Attached_list2(i)==Charged_ID
                   total2 = total2-1;   % Protonation site
                   total4 = total4-1;
                   Charge2=Charge2+1;
                   Charge4=Charge4+1;
                end
            end    
            Mul1 = rem(total1,2)+1; Mul2 = rem(total2,2)+1;
            Mul3 = rem(total3,2)+1; Mul4 = rem(total4,2)+1;
            WritePartCom(directory, [filename1,'.com'], Atom_info, Attached_list1, Charge1, Mul1);
            WritePartCom(directory, [filename2,'.com'], Atom_info, Attached_list2, Charge2, Mul2);
            WritePartCom(directory, [filename3,'.com'], Atom_info, Attached_list1, Charge3, Mul3);
            WritePartCom(directory, [filename4,'.com'], Atom_info, Attached_list2, Charge4, Mul4);
    end
    
    function WritePartCom(directory, filename, Atom_info, Attached_list, Charge, Mul)
        fileID = fopen(fullfile(directory,filename),'w');
        fprintf(fileID,'%%NProcShared=12\n');
        fprintf(fileID,['%%Chk=',filename,'.chk\n']);
        fprintf(fileID,'#n B3LYP/6-31G(d) Opt Freq pressure=0.00001\n\n');
        fprintf(fileID, [filename,'\n\n']);
        fprintf(fileID, [num2str(Charge),' ',num2str(Mul),'\n']);
        for i = 1:length(Attached_list)
            index = Attached_list(i);
            fprintf(fileID,[Atom_info(index).element,'\t']);
            fprintf(fileID,'%.5f\t',Atom_info(index).X);
            fprintf(fileID,'%.5f\t',Atom_info(index).Y);
            fprintf(fileID,'%.5f\n',Atom_info(index).Z);
        end
        fprintf(fileID,'\n');
        fclose(fileID);
    end
end