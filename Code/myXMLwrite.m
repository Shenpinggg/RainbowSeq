function ST = myXMLwrite(fileName, docNode)
% ST = myXMLwrite(FILENAME, docNode) Saves (properly!) an XML file.
%   default command XMLWRITE causes addition of extra line breaks in the xml
%   file, this version worksaround this apparent bug
%   FILENAME is the file name string to save to (e.g. 'test.xml')
%   DOCNODE is a Document Object Model (most commonly created by the
%   XMLREAD command, this object will be written to the file specified by
%   FILENAME.
%
%   The function returns success status from the FCLOSE function - returns
%   0 if successful, -1 if not.
%
%		EXAMPLE:
%			docNode=xmlread('somefile.xml');
%			% Do XML stuff here... %
%			myXMLwrite('somefile.xml', docNode);
%
%
%		Created by Gilad Yahalom
%		Fell free to modify code, please leave acknowledgments to original author
%
%   See also xmlread, xmlwrite.
%

docNodeRoot=docNode.getDocumentElement;
str=strrep(char(docNode.saveXML(docNodeRoot)), 'encoding="UTF-16"', 'encoding="UTF-8"');
fid=fopen(fileName, 'w');
fwrite(fid, str);
ST = fclose(fid);
