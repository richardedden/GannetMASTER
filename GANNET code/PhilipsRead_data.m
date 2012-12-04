function [ MRS_struct ] = PhilipsRead_data(MRS_struct, fname, fname_water )
% RE/CJE Parse SPAR file for header info
% 110825
   
   % work out data header name
   sparname = [fname(1:(end-4)) 'list'];
   sparheader = textread(sparname, '%s');
   sparidx=find(ismember(sparheader, 'F-resolution')==1);
   MRS_struct.npoints = str2num(sparheader{sparidx+2});
   sparidx=find(ismember(sparheader, 'number_of_extra_attribute_1_values')==1);
   MRS_struct.nrows = str2num(sparheader{sparidx+2});
   
   sparidx=find(ismember(sparheader, 'number_of_signal_averages')==1);
   %MRS_struct.Navg = MRS_struct.nrows * str2num(sparheader{sparidx+2});
   MRS_struct.Navg(MRS_struct.ii) = str2num(sparheader{sparidx+2}); %Trial SDAT might be average not sum.
   %sparidx=find(ismember(sparheader, 'repetition_time')==1);
   %MRS_struct.TR = MRS_struct.nrows * str2num(sparheader{sparidx+2});
   %TR not contained in .data - hard code for now.
   MRS_struct.TR = 2.0;
   
   %sparidx=find(ismember(sparheader, 'sample_frequency')==1);
   %MRS_struct.sw = str2num(sparheader{sparidx+2});
    %SW not contained in .data - hard code for now.
   MRS_struct.sw = 2000;
   
   
   %Need to determine the offset of the data - i.e. how many channels are
   %there...
   sparidx=find(ismember(sparheader, 'NOI')==1);
   MRS_struct.coil_chammels=size(sparidx,1)-2;
   sparidx=find(ismember(sparheader, 'STD')==1);
   MRS_struct.ptr_offset=str2num(sparheader{sparidx(3)+20});
   
   
   %MRS_struct.data = SDATreadMEGA(fname, MRS_struct.npoints, MRS_struct.nrows);
   %Need to skip rows associated with the '
   MRS_struct.data = readraw_Gannet(fname, 'float', [2 MRS_struct.npoints 1 MRS_struct.Navg(MRS_struct.ii) MRS_struct.nrows], 'l',MRS_struct.ptr_offset);
   
   %  Make data complex.
   MRS_struct.data = squeeze(MRS_struct.data(1,:,:,:,:)+ 1i*MRS_struct.data(2,:,:,:,:));
   
   
   
   %undo time series phase cycling to match GE
   corrph = ones(size(MRS_struct.data));
   for jj=1:size(MRS_struct.data,3)
    corrph(:,:,jj) = corrph(:,:,jj) * (-1).^(jj+1);
   end
   
   MRS_struct.data = MRS_struct.data .* corrph;
   %MRS_struct.data = MRS_struct.data .*repmat(conj(MRS_struct.data(1,:))./abs(MRS_struct.data(1,:)),[MRS_struct.npoints 1]);
   %Philips data appear to be phased already (ideal case)
   MRS_struct.data = conj(MRS_struct.data); %RE 110728 - empirical factor to scale 'like GE'
   MRS_struct.data = reshape(MRS_struct.data,[size(MRS_struct.data,1) size(MRS_struct.data,2)*size(MRS_struct.data,3)]);
   % Depending on the ordering of OFF and ON, the minus sign here may have
   % to be removed.... not sure how to automate/fix... raee 4/9/12
   MRS_struct.data = MRS_struct.data .*repmat(conj(MRS_struct.data(1,:))./abs(MRS_struct.data(1,:)),[MRS_struct.npoints 1]);
end

