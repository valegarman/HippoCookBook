function[ExtractedPhotometryAndBehaviour]=ExtractPhotometryWithoutBehaviour(PhotometryFileName)

    data_struct = import_ppd(PhotometryFileName);
    [dFF,Signal1_denoised,Signal2_denoised,Signal1_normalized,Signal2_normalized] = ComputeDFF(data_struct);
    TriggeredEvents=find(diff(data_struct.digital_1)==1);
    CameraFrameEvents=find(diff(data_struct.digital_2)==1);
   
   BehaviourStartEndCropping =[0 0];
    dFFFiltered=nanmean([dFF(CameraFrameEvents),dFF(max(1,CameraFrameEvents-1)),dFF(CameraFrameEvents+1)],2);
    dFFFilteredCropped=dFFFiltered(1-BehaviourStartEndCropping(1):end-BehaviourStartEndCropping(2));
%          TriggeredFiltered=(TriggeredEvents- CameraFrameEvents(1,1));

      TriggeredFilteredp=data_struct.digital_1(CameraFrameEvents);
      TriggeredFiltered=find(diff(TriggeredFilteredp)==1);

    ZScoreFiltered=(dFFFilteredCropped-nanmean(dFFFilteredCropped))/nanstd(dFFFilteredCropped);
    ZScoreCropped=ZScoreFiltered;
    
    ExtractedPhotometryAndBehaviour.Photometry.dFF=dFF;
    ExtractedPhotometryAndBehaviour.Photometry.Signal1_normalized=Signal1_normalized;
    ExtractedPhotometryAndBehaviour.Photometry.Signal2_normalized=Signal2_normalized;
    ExtractedPhotometryAndBehaviour.Photometry.Signal1_denoised=Signal1_denoised;
    ExtractedPhotometryAndBehaviour.Photometry.Signal2_denoised=Signal2_denoised;

    ExtractedPhotometryAndBehaviour.Photometry.TriggeredEvents=TriggeredEvents;
    ExtractedPhotometryAndBehaviour.Photometry.CameraFrameEvents=CameraFrameEvents;
    ExtractedPhotometryAndBehaviour.Photometry.TriggeredFiltered=TriggeredFiltered;

       ExtractedPhotometryAndBehaviour.Behaviour.dFFFilteredCropped=dFFFilteredCropped;
    ExtractedPhotometryAndBehaviour.Behaviour.ZScoreCropped=ZScoreCropped;
    
end