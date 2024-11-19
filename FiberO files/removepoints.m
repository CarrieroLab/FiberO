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
   
   This part of the FiberO software allows the manual removal of selected
   regions by the user

   Inputs:
       img          - Grayscale image.
       removedback  - Tissue image after background removal.
       average      - Average intensity of each "x_y_region".

   Outputs:
       img          - The updated image after the user has removed selected regions.
       removedback  - Tissue image after background removal.
       Selection    - User-selected regions for specific area removal.

%}
function [img,removedback,Selection]=removepoints(img,removedback,average)

figure()
imagesc(img); colormap gray; hold on; axis off; title('Removal selection'); truesize;

ngrid=100;
g = size(img,1)/ngrid;
Selection=ones(size(img));
key = waitforbuttonpress;

while key == 0
roi=roipoly;
img=double(img).*double(~roi);
Selection=Selection.*double(~roi);
imagesc(img); colormap gray;
key = waitforbuttonpress;
end

close(gcf)

end