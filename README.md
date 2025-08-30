# Wirewalker
Nortek Signature 1000 and RBR Data processing for the Wirewalker vehicle.
Maintained by MOD/SIO for Scripps Users. This processing will take some time to get set up for new users.

Processing of the WireWalker data using the upcast data:
- processing of the aquadopp data - from beam velocity to geographical
velocities (North-South-Up)

- Read ctd data (Ruskin). Define upcasts grid aquadopp and ctd data on
the same vertical-time space

- Attempt to correct/filter the surface wave filter. 2 options :
wavelet filter of the upcasts or 2d filter of the gridded product.

Watch out: for some experiment the RBR ctd did not store the data correctly. Channels are mess up but the data is still good. 
In that case contact us or check out Toolbox/RSK_structM.m

It is very important to check out the path to the data since it also change with the experiment. 
This "meta-data" aspect of the processing will be subject to improvement in the futur.  

Feed back are welcome to help us improvements


Developed by D.Lucas and A.Le Boyer
10/28/2016
