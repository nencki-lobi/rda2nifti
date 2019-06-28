
%
% Read spectroscopy data from Siemens machine
%
% Read a .rda file
%
%
if ~exist('rda_filename','var')
[filename , pathname ] = uigetfile('*.rda', 'Select an RDA file')
rda_filename = [pathname , filename]; %'c:/data/spectroscopy/spec raw data/MrSpec.20020531.160701.rda'
end
fid = fopen(rda_filename);

head_start_text = '>>> Begin of header <<<';
head_end_text   = '>>> End of header <<<';

tline = fgets(fid)

while (isempty(strfind(tline , head_end_text)))
    
    tline = fgets(fid)
    
    if ( isempty(strfind (tline , head_start_text)) + isempty(strfind (tline , head_end_text )) == 2)
        
        
        % Store this data in the appropriate format
        
        occurence_of_colon = findstr(':',tline);
        variable = tline(1:occurence_of_colon-1) ;
        value    = tline(occurence_of_colon+1 : length(tline)) ;
        
        switch variable
        case { 'PatientID' , 'PatientName' , 'StudyDescription' , 'PatientBirthDate' , 'StudyDate' , 'StudyTime' , 'PatientAge' , 'SeriesDate' , ...
                    'SeriesTime' , 'SeriesDescription' , 'ProtocolName' , 'PatientPosition' , 'ModelName' , 'StationName' , 'InstitutionName' , ...
                    'DeviceSerialNumber', 'InstanceDate' , 'InstanceTime' , 'InstanceComments' , 'SequenceName' , 'SequenceDescription' , 'Nucleus' ,...
                    'TransmitCoil' }
            eval(['rda.' , variable , ' = value ']);
        case { 'PatientSex' }
            % Sex converter! (int to M,F,U)
            switch value
            case 0
                rda.sex = 'Unknown';
            case 1
                rda.sex = 'Male';
            case 2
                
                rda.sex = 'Female';
            end
            
        case {  'SeriesNumber' , 'InstanceNumber' , 'AcquisitionNumber' , 'NumOfPhaseEncodingSteps' , 'NumberOfRows' , 'NumberOfColumns' , 'VectorSize' }
            %Integers
            eval(['rda.' , variable , ' = str2num(value) ']);
        case { 'PatientWeight' , 'TR' , 'TE' , 'TM' , 'DwellTime' , 'NumberOfAverages' , 'MRFrequency' , 'MagneticFieldStrength' , 'FlipAngle' , ...
                     'SliceThickness' ,  'FoVHeight' , 'FoVWidth' , 'PercentOfRectFoV' , 'PixelSpacingRow' , 'PixelSpacingCol'}
            %Floats 
            eval(['rda.' , variable , ' = str2num(value) ']);
        case {'SoftwareVersion[0]' }
            rda.software_version = value;
        case {'CSIMatrixSize[0]' }
            rda.CSIMatrix_Size(1) = str2num(value);    
        case {'CSIMatrixSize[1]' }
            rda.CSIMatrix_Size(2) = str2num(value);    
        case {'CSIMatrixSize[2]' }
            rda.CSIMatrix_Size(3) = str2num(value);    
        case {'PositionVector[0]' }
            rda.PositionVector(1) = str2num(value);    
        case {'PositionVector[1]' }
            rda.PositionVector(2) = str2num(value);     
        case {'PositionVector[2]' }
            rda.PositionVector(3) = str2num(value);    
        case {'RowVector[0]' }
            rda.RowVector(1) = str2num(value);    
        case {'RowVector[1]' }
            rda.RowVector(2) = str2num(value);       
        case {'RowVector[2]' }
            rda.RowVector(3) = str2num(value);    
        case {'ColumnVector[0]' }
            rda.ColumnVector(1) = str2num(value);     
        case {'ColumnVector[1]' }
            rda.ColumnVector(2) = str2num(value);       
        case {'ColumnVector[2]' }
            rda.ColumnVector(3) = str2num(value);  
                    case {'VOINormalCor' }
            rda.VOINormalCor = str2num(value);  
                    case {'VOINormalSag' }
            rda.VOINormalSag = str2num(value);  
                    case {'VOINormalTra' }
            rda.VOINormalTra = str2num(value);  
                    case {'VOIPositionCor' }
            rda.VOIPositionCor = str2num(value);  
                    case {'VOIPositionSag' }
            rda.VOIPositionSag = str2num(value);  
                    case {'VOIPositionTra' }
            rda.VOIPositionTra = str2num(value);  
                    case {'VOIRotationInPlane' }
            rda.VOIRotationInPlane = str2num(value);  
            %VOIRotationInPlane

            
        otherwise
            % We don't know what this variable is.  Report this just to keep things clear
            disp(['Unrecognised variable ' , variable ]);
        end
        
    else
        % Don't bother storing this bit of the output
    end
    
end

%
% So now we should have got to the point after the header text
% 
% Siemens documentation suggests that the data should be in a double complex format (8bytes for real, and 8 for imaginary?)
%

bytes_per_point = 16;
complex_data = fread(fid , rda.CSIMatrix_Size(1) * rda.CSIMatrix_Size(1) *rda.CSIMatrix_Size(1) *rda.VectorSize * 2 , 'double');  

%fread(fid , 1, 'double');  %This was a check to confirm that we had read all the data (it passed!)
fclose(fid);

% Now convert this data into something meaningful

 %Reshape so that we can get the real and imaginary separated
 hmm = reshape(complex_data,  2 , rda.VectorSize , rda.CSIMatrix_Size(1) ,  rda.CSIMatrix_Size(2) ,  rda.CSIMatrix_Size(3) );
 
 %Combine the real and imaginary into the complex matrix
 hmm_complex = complex(hmm(1,:,:,:,:),hmm(2,:,:,:,:));
 
 %Remove the redundant first element in the array
% % %  Time_domain_data = reshape(hmm_complex, rda.VectorSize , rda.CSIMatrix_Size(1) ,  rda.CSIMatrix_Size(2) ,  rda.CSIMatrix_Size(3));
% % % 
% % %  
% % %  Time_domain_data = fft(fft(fftshift(fftshift(padarray(fftshift(fftshift(ifft(ifft(Time_domain_data,[],2),[],3),2),3),[0 16 16],'both'),2),3),[],2),[],3);
% % %       
% % % % Time_domain_data = ifft(fftshift(fft(ifft(fftshift(fft(Time_domain_data,[],2),2),64,2),[],3),3),64,3);
% % %  
% % %  
% % %  %Time_domain_data = ifft(fftshift(fft(ifft(fftshift(fft(Time_domain_data,[],2),2),64,2),[],3),3),64,3);
% % % 
% % %  %Time_domain_data = fft(fft(fftshift(fftshift(ifft(ifft(Time_domain_data,[],2),[],3),2),3),64,2),64,3);
% % % %
% % % % Do a crude image
% % % for counter = 1 : 100
% % %     figure(4)
% % %     T2 = 200;
% % %     weighting_function = exp(-[10:10:640]'/T2);
% % %     weighted_time_domain = Time_domain_data;
% % %     for x_counter = 1:size(Time_domain_data,2)
% % %         for y_counter = 1: size(Time_domain_data,3)
% % %             weighted_time_domain(:,x_counter,y_counter) = weighting_function .* Time_domain_data(:,x_counter,y_counter);
% % %         end
% % %     end
% % %     
% % %  imagesc(interp2(interp2(squeeze(abs(mean(weighted_time_domain(:,:,:)))))'));
% % %  axis('ij')
% % % figure(1)
% % %  imagesc(squeeze(angle(mean(Time_domain_data(2:5,:,:))))');
% % %  axis('ij')
% % %  
% % %  %Get a pixel and show a spectrum on that point
% % %  disp('Select pixel for a spectrum');
% % % [xtmp,ytmp] = getpts
% % % 
% % % if (length(xtmp) >= 2)
% % %     %create box with the range of min(xtmp) to max(xtmp) and min(ytmp) to
% % %     %max(ytmp)
% % %     
% % %     flag = 0;
% % %     for x_pos = floor(min(xtmp+0.5)) : floor(max(xtmp+0.5))
% % %         for y_pos = floor(min(ytmp+0.5)) : floor(max(ytmp+0.5))
% % %             if (flag == 0)
% % %             time_domain = (fftshift(fft(Time_domain_data(:,floor(x_pos),floor(y_pos)))));
% % %             fid=  Time_domain_data(:,floor(x_pos),floor(y_pos));
% % %             flag = 1;
% % %         else
% % %             time_domain = time_domain + (fftshift(fft(Time_domain_data(:,floor(x_pos),floor(y_pos)))));
% % %             fid= fid + Time_domain_data(:,floor(x_pos),floor(y_pos));
% % %         end
% % %         end
% % %     end
% % %      figure(2)
% % %  %show the spectrum of that pixel or region
% % %  
% % %  plot (abs(time_domain));
% % %  
% % %   figure(3)
% % %  %show the fid of that pixel
% % %  
% % %  plot (real(fid));
% % %  hold on
% % %  plot(imag(fid),'r-');
% % %  hold off
% % % 
% % % else
% % %  figure(2)
% % %  %show the spectrum of that pixel or region
% % %  
% % % plot (abs(fftshift(fft(Time_domain_data(:,floor(xtmp(1)),floor(ytmp(1)))))));
% % %  
% % %   figure(3)
% % %  %show the fid of that pixel
% % %  
% % %  plot (real(((Time_domain_data(:,floor(xtmp(1)),floor(ytmp(1)))))));
% % %  hold on
% % %  plot (imag(((Time_domain_data(:,floor(xtmp(1)),floor(ytmp(1)))))),'r-');
% % %  hold off
% % % end
% % % end