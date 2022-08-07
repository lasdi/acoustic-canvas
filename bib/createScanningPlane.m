function [thetaScanAngles,phiScanAngles,scanningPointsX,...
    scanningPointsY,scanningPointsZ] = createScanningPlane(...
    maxScanningPlaneExtentL,maxScanningPlaneExtentH,numberOfScanningPointsL,...
    numberOfScanningPointsH,distanceToScanningPlane,xPos,yPos,zPos)
%CREATESCANNINGPLANE Summary of this function goes here
%   Detailed explanation goes here
% Create scanning points
scanningAxisX = -maxScanningPlaneExtentL:2*maxScanningPlaneExtentL/(numberOfScanningPointsL-1):maxScanningPlaneExtentL;
scanningAxisY = -maxScanningPlaneExtentH:2*maxScanningPlaneExtentH/(numberOfScanningPointsH-1):maxScanningPlaneExtentH;
[scanningPointsX, scanningPointsY] = meshgrid(scanningAxisX, scanningAxisY);



%Get angles to scanning points and source positions
[thetaScanAngles, phiScanAngles] = convertCartesianToSpherical(scanningPointsX, scanningPointsY, distanceToScanningPlane);

% Adjust for plot
scanningPointsZ =ones(size(scanningPointsX))*distanceToScanningPlane;
scanningPointsX = scanningPointsX + mean(xPos);
scanningPointsY = scanningPointsY + mean(yPos);
scanningPointsZ = scanningPointsZ + mean(zPos);
end

