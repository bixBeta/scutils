for i in batch_*
do
	cd $i
	cp ../sobj_report.qmd .

	quarto render sobj_report.qmd -P title:${i}_Report -o ${i}_Report.html

	rm sobj_report.qmd

	cd ..

done
