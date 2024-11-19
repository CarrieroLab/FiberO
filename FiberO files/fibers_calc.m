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
   
   This part of the FiberO software specifically deals with calculating fiber
   organization across the tissue.

   Inputs:
       img        - Grayscale image.
       h          - Contains the orientation value (in deg) of each "x_y_region".
       Selection  - User-selected regions for specific area removal.
       Tcr        - Threshold for fibers.
   
   Outputs:
       B               - Mask containing the tissue surfaces with continuous fibers.
       C               - Mask containing the fibers continuity.
       cont            - Index for posterior classification of surfaces.
       indexOfBiggest  - Numbers of groups of fibers surfaces.

%}
function [B, C, cont, indexOfBiggest]=fibers_calc(img, h, Selection, Tcr)

for i = 2:size(h,1)-1
    for j = 2:size(h,1)-1
        l = zeros(size(h));
        l(i,j) = 1;
        neighbors = h(conv2(l,[1,1,1;1,0,1;1,1,1],'same')>0);
        df = h(i,j) - neighbors;
        cnnctn = find(abs(df)<20 | abs(df)>=160);                          
        n{i,j} = cnnctn;
    end
end

n(size(h,1),size(h,2)) = {[]};
M = zeros(size(n));
chr = 0;
cont = 0;
arry = zeros(1,8);

for i = 1:size(n,1)
    for j = 1:size(n,2)
        pos = [i-1, j-1; i, j-1; i+1, j-1; i-1, j; i+1, j; i-1, j+1; i, j+1; i+1, j+1];     
        nmbrs = n{i,j};                                                                     
        for p = 1:length(nmbrs)
            for u = 1:length(nmbrs)
                if M(pos(nmbrs(u),1),pos(nmbrs(u),2)) == 0
                else
                    arry(u) = M(pos(nmbrs(u),1), pos(nmbrs(u),2));
                end
            end
            arry = [arry M(i,j)];                                          
            chr = min((arry(find(arry~=0)))); 
            chr2 = max((arry(find(arry~=0)))); 
            if isempty(chr); chr = 0; chr2 = 0;end
            if M(i,j) ~= 0
                if chr ~= chr2
                    M(M == chr2) = chr;
                end
            else
                if chr ~= chr2
                    M(i,j) = chr;
                    M(M == chr2) = chr;
                else
                    cont = cont + 1;
                    M(i,j) = cont;
                end
            end
            if chr == 0
                M(pos(nmbrs(p),1), pos(nmbrs(p),2)) = cont;
            else
                M(pos(nmbrs(p),1), pos(nmbrs(p),2)) = chr;                 
            end
        end
        arry = zeros(1,8);
    end
end

allAreas = regionprops(M,'Area');
allAreas=[allAreas.Area];
[biggestArea, indexOfBiggest] = sort(allAreas,'descend');                  
indexOfBiggest = indexOfBiggest(find(biggestArea(biggestArea>=30)));      
C = zeros(size(img)); B = C;
cont = 2.1;

for i = 1:size(indexOfBiggest,2)
    blob = M == indexOfBiggest(i);                                         
    blob = M.*blob;
    blob=padarray(blob,[1 1], 'both');                                     
    for j=2:size(blob,1)-1
        for k=2:size(blob,1)-1
            pos = [j-1, k-1; j, k-1; j+1, k-1; j-1, k; j+1, k; j-1, k+1; j, k+1; j+1, k+1];
            T=0;
            if blob(j,k)== 0
                for v=1:size(pos,1)
                    T=sum(blob(pos(v,1),pos(v,2)))+T;
                end
                if T>indexOfBiggest(i)*4                                   
                    blob(j,k)=indexOfBiggest(i);
                end
            end
        end
    end
    blob=blob(2:size(blob,1)-1,2:size(blob,2)-1);                          
    conn = imfill(imresize(blob,size(img)/size(M),'box'));                 
    conn = conn.*Selection;
    conns{1,i} = conn;
    if size(conns,2)>1
        for j = 1:size(conns,2)-1
            idx1 = [];
            idx2 = [];
            AA = logical(conn);
            BB = logical(conns{1,j});
            for a=1:size(conn,1)
                for b=1:size(conn,2)
                    if AA(a,b)==BB(a,b) && AA(a,b)~=0 && BB(a,b) ~=0
                        conn(a,b) = 0;
                    end
                end
            end
        end
    end
    conns{i} = conn;
    dati = double(img).*conn;
    srfc_area = sum(sum(conn~=0));
    conn(conn~=0) = cont;
    cont = cont + 1;
    C = C + conn;
    sizesrfc(i) = srfc_area/sum(sum(img>Tcr));
    maskang=h.*blob>0;
    maskang=h(maskang>0);
    for j=1:length(maskang)
        if maskang(j)<90
            maskang(j)=maskang(j)+180;
        end
    end
    if sizesrfc(i)<1/6; continue; end
    B = B + conn;
end
B(B>1)=1;
end