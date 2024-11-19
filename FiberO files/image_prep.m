%{
   Carriero Lab - City College of New York

   FiberO to Quantify Fibers Orientation and Organization
   FiberO is a software designed to quantify fibers orientation by performing 
   morphological openings and automatically measure fibers organization by 
   determining and plotting fiber continuity across the tissue

   Conceptualized by: Alessandra Carriero
   Written by: Asier Munoz, Anxhela Docaj, Julen Fernandez

   LICENSE:
   FiberO is a free software Licensed under the GNU Affero General Public License v3.0 License 
   as published by the Free Software Foundation.
   You may not use this file except in compliance with the License.
   You may obtain a copy of the License in the FiberO folder or at

   https://www.gnu.org/licenses/agpl-3.0.en.html 

   Unless required by applicable law or agreed to in writing, software
   distributed under this License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   This part of the FiberO software specifically deals with image
   preprocessing for fiber enhancement

   Inputs:
       img1  - Image used for preprocessing.
       Tr    - Threshold for tissue.
   
   Outputs:
       img          - Preprocessed image used further steps.
       imgo         - Image used for the calculation of tissue area.
       x_y_region   - Created grid.
       average      - Average intensity of each "x_y_region".
       removedback  - Tissue image after background removal.
       img0         - Original unprocessed image.

%}
function [img, imgo, x_y_region, average, removedback, img0] = image_prep(img1,Tr)

if size(img1,1)==512
    img1 = img1(5:size(img1,1)-8,7:size(img1,2)-6,:);
end
if size(img1,1)==1024
    img1 = img1(9:size(img1,1)-16,13:size(img1,2)-12,:);
end
if size(img1,1)==2048
    img1 = img1(17:size(img1,1)-32,25:size(img1,2)-24,:);
end

% For a resolution different to 512,1024 or 2048
if size(img1,1)>size(img1,2)
    img1 = img1(round((size(img1,1)-size(img1,2))/2):round(size(img1,2)+(size(img1,1)-size(img1,2))/2)-1,1:size(img1,2));
end
if size(img1,1)<size(img1,2)
    img1 = img1(1:size(img1,1),round((size(img1,2)-size(img1,1))/2):round(size(img1,1)+(size(img1,2)-size(img1,1))/2)-1);
end

imgo=img1;
img0=img1;

I2_1=imfilter(img1,(1/9)*ones(3,3),'conv');

Gx=imfilter(I2_1,[-1 -2 -1;0 0 0 ;1 2 1],'conv');
Gy=imfilter(I2_1,[-1 0 1;-2 0 2 ;-1 0 1],'conv');
img1=uint16(I2_1+((Gx+Gy)/2));


ngrid = 100;                                                    
g = size(img1,2)/ngrid;                                       
hstimng = img1;
fibermask=imbinarize(hstimng, 'adaptive', 'Sensitivity', 0.95, 'ForegroundPolarity', 'bright');
img=double(hstimng).*fibermask;

posGrid = 0;
for i = 1:ngrid                                                  
    if (posGrid < size(hstimng,1))
        x_y_region(i,:) = [g*(i-1) g*i];
    else
        break
    end
    posGrid = g*i;
end
x_y_region(1,1) = x_y_region(1,1) + 1;
x_y_region = floor(x_y_region);

for i=1:size(x_y_region,1)
    for j=1:size(x_y_region,1)
        average(i,j)=mean(mean(imgo(x_y_region(i,1):x_y_region(i,2),x_y_region(j,1):x_y_region(j,2))));
        if average(i,j)<Tr         
            imgo(x_y_region(i,1):x_y_region(i,2),x_y_region(j,1):x_y_region(j,2))=0;
            average(i,j)=0;
        end
    end
end
removedback=imgo;
for i=1:size(x_y_region,1)
    for j=1:size(x_y_region,1)
        average(i,j)=mean(mean(img(x_y_region(i,1):x_y_region(i,2),x_y_region(j,1):x_y_region(j,2))));
        if average(i,j)<Tr
            img(x_y_region(i,1):x_y_region(i,2),x_y_region(j,1):x_y_region(j,2))=0;
            average(i,j)=0;
        else
            average(i,j)=1;
        end
    end
end
end