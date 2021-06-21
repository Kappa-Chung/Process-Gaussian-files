function [Atom, Bond, Energy] = Read_Log(directory, filename)
    % Input directory and filename of gaussian .log file
    % The function will output the information of atoms, bond, and energies
    % written in the log file.

    file = fullfile(directory,filename);
    fileID = fopen(file, 'r');  % Open the file, with read permission
    exp1 = '[^\s:(),]*'; % Expression without space and colon

    keySet = {'1','6','7','8','9','15','16','17'};
    valueSet = {'H','C','N','O','F','P','S','Cl'};
    AtomMap = containers.Map(keySet,valueSet);
    
    % Create empty struct for bond
    fields = {'index','Atom1','Atom2','length','Wiberg'}; 
    Bond = cell2struct(cell(5,1),fields);
    
    
    while feof(fileID)==0
        line = fgetl(fileID);
        b1=regexp(line,exp1,'match');
        if size(b1,2)>=2
            switch [b1{1},' ',b1{2}]
                case {'Standard orientation'}
                    for i = 1:4
                        line = fgetl(fileID);
                    end
                    i=1;
                    while feof(fileID)==0
                        line = fgetl(fileID);
                        b1=regexp(line,exp1,'match');
                        if strcmp(b1{1}(1),'-')
                           break 
                        end
                        Atom(i).index = str2double(b1{1});
                        Atom(i).element = AtomMap(b1{2});
                        Atom(i).type = str2double(b1{3});
                        Atom(i).X = str2double(b1{4});
                        Atom(i).Y = str2double(b1{5});
                        Atom(i).Z = str2double(b1{6});
                        i=i+1;
                    end
                case {'Mulliken charges'}
                    line = fgetl(fileID);
                    for i = 1:size(Atom,2)
                        line = fgetl(fileID);
                        b1=regexp(line,exp1,'match');
                        index = find([Atom.index]==str2double(b1{1}));
                        Atom(index).Mulliken = str2double(b1{3});
                    end
                    for i = 1:3
                        line = fgetl(fileID);
                    end
                    for i = 1:sum(~strcmp({Atom.element},'H'))
                        line = fgetl(fileID);
                        b1=regexp(line,exp1,'match');
                        index = find([Atom.index]==str2num(b1{1}));
                        Atom(index).Mulliken_Hsum = str2double(b1{3});
                    end
                case {'APT charges'}
                    line = fgetl(fileID);
                    for i = 1:size(Atom,2)
                        line = fgetl(fileID);
                        b1=regexp(line,exp1,'match');
                        index = find([Atom.index]==str2double(b1{1}));
                        Atom(index).APT = str2double(b1{3});
                    end
                    for i = 1:3
                        line = fgetl(fileID);
                    end
                    for i = 1:sum(~strcmp({Atom.element},'H'))
                        line = fgetl(fileID);
                        b1=regexp(line,exp1,'match');
                        index = find([Atom.index]==str2num(b1{1}));
                        Atom(index).APT_Hsum = str2double(b1{3});
                    end
                case {'Wiberg bond'}
                    for j = 1:ceil(size(Atom,2)/9)
                        for i = 1:3
                            line = fgetl(fileID);
                        end

                        for i = 1:size(Atom,2)
                            line = fgetl(fileID);
                            b1=regexp(line,exp1,'match');
                            for k = 3:length(b1)
                                index = [[Bond.Atom1]==i & [Bond.Atom2]==((j-1)*9)+k-2];
                                if any(index)
                                    Bond(index).Wiberg = str2num(b1{k});
                                end
                            end
                        end 
                    end
                    
                    for i = 1:6
                        line = fgetl(fileID);
                    end
                    for i = 1:size(Atom,2)
                        line = fgetl(fileID);
                        b1=regexp(line,exp1,'match');
                        index = find([Atom.index]==i);
                        Atom(index).Wiberg = str2double(b1{3});
                    end
                    
                case {'Thermal correction'}
                    for i = 1:3
                        line = fgetl(fileID);
                    end
                    b1=regexp(line,exp1,'match');
                    Energy.Electronic_ZeroPoint_Energy = str2double(b1{end});
                    line = fgetl(fileID);
                    b1=regexp(line,exp1,'match');
                    Energy.Electronic_Thermo_Energy = str2double(b1{end});
                    line = fgetl(fileID);
                    b1=regexp(line,exp1,'match');
                    Energy.Electronic_Thermo_Enthalpy = str2double(b1{end});
                    line = fgetl(fileID);
                    b1=regexp(line,exp1,'match');
                    Energy.Electronic_ThermoFree_Energy = str2double(b1{end});
                case {'! Optimized'}
                    for i = 1:4
                        line = fgetl(fileID);
                    end
                    i=1;
                    while feof(fileID)==0
                        line = fgetl(fileID);
                        b1=regexp(line,exp1,'match');
                        if strcmp(b1{1}(1),'-') % In case only R (two atoms)
                           break 
                        end
                        if ~strcmp(b1{2}(1),'R')
                           break 
                        end
                        Bond(i).index = str2double(b1{2}(2:end));
                        Bond(i).Atom1 = str2double(b1{4});
                        Bond(i).Atom2 = str2double(b1{5});
                        Bond(i).length = str2double(b1{6});
                        i=i+1;
                    end
            end
        end
        
        if ~isempty(b1)
            if any(strcmp(b1{end}(end),'@'))
                    line = fgetl(fileID);
                    line = fgetl(fileID);            
                    while feof(fileID)==0
                        line = fgetl(fileID);
                        b1=regexp(line,exp1,'match');
                        if ~isempty(b1)
                            if strcmp(b1{1},'Job')
                                disp(' ')
                                break
                            end
                        end
                        disp(line)
                    end
            end
        end
    end
    fclose(fileID); 
end