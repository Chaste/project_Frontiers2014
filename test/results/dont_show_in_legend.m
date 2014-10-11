function dont_show_in_legend(object_handle, true_or_false)
hAnnotation = get(object_handle,'Annotation');
hLegendEntry = get(hAnnotation','LegendInformation');
if (true_or_false)
    set(hLegendEntry,'IconDisplayStyle','off')
end